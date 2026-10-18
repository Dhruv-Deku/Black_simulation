#include "bodies.h"

#include "blackbody.h"
#include "central.h"
#include "gravity.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {

constexpr double kHalfPi = 0.5 * kerr::kPi;
constexpr size_t kMaxTrail = 3000;
constexpr size_t kMaxEvents = 8;
constexpr size_t kMaxImpacts = 16;
constexpr int kBackfillSamples = 1400;  // 700 M of past orbit for each new body (times the scene scale)
// A tidally stretched star is treated as torn apart once its longest axis reaches
// this many radii: from there on its pieces are effectively unbound.
constexpr double kDisruptionStretch = 2.5;
// Gravity between bodies is softened over this fraction of their radii.
constexpr double kSoftening = 0.1;
// A small black hole bends light out to this many times its mass; beyond, the bend
// (4m/b < 0.002 rad) is below a pixel and is faded out.
constexpr double kLensReach = 2000.0;

std::string format(const char* fmt, double a, double b = 0.0) {
    char buf[200];
    std::snprintf(buf, sizeof(buf), fmt, a, b);
    return buf;
}

// Position (x, y, z, r), 4-velocity and shape, as stored in the history.
History::Sample sampleOf(double spin, const Body& b) {
    const kerr::State& s = b.state;
    auto p = kerr::toCartesian(spin, s.r, s.theta, s.phi);
    auto v = kerr::fourVelocity(spin, b.constants, s);
    History::Sample out;
    out.position = glm::vec4(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]),
                             static_cast<float>(s.r));
    out.u = glm::vec4(static_cast<float>(v[0]), static_cast<float>(v[1]), static_cast<float>(v[2]),
                      static_cast<float>(v[3]));
    if (b.kind == BodyKind::TidalStar) {
        const tidal::Mat3 S = tidal::inverseShape(b.shape.F);
        out.shapeA = glm::vec4(static_cast<float>(S[0][0]), static_cast<float>(S[1][1]), static_cast<float>(S[2][2]),
                               static_cast<float>(tidal::stretch(b.shape.F)));
        out.shapeB = glm::vec4(static_cast<float>(S[0][1]), static_cast<float>(S[0][2]), static_cast<float>(S[1][2]), 0.0f);
    } else {
        out.shapeA = glm::vec4(1.0f);
        out.shapeB = glm::vec4(0.0f);
    }
    return out;
}

History::Sample mix(const History::Sample& a, const History::Sample& b, float f) {
    return {glm::mix(a.position, b.position, f), glm::mix(a.u, b.u, f), glm::mix(a.shapeA, b.shapeA, f),
            glm::mix(a.shapeB, b.shapeB, f)};
}

int column(long long k) {
    return static_cast<int>(((k % History::kColumns) + History::kColumns) % History::kColumns);
}

// d(tau)/dt along a geodesic: Sigma / (dt/dMino).
double clockRate(double spin, const kerr::Constants& c, const kerr::State& s) {
    return kerr::metricAt(spin, s.r, s.theta).sigma / kerr::derivatives(spin, c, s).dt;
}

gravity::Member memberOf(const Body& b) {
    return {b.constants, b.state, b.mass, kSoftening * b.radius};
}

// Lorentz factor of a body relative to the rigidly rotating surface it strikes:
// -u_body . u_surface, with u_surface = u^t (1, 0, 0, Omega).
double impactGamma(double spin, const Body& b, double omega) {
    const kerr::Vec4 u = kerr::fourVelocity(spin, b.constants, b.state);
    const double ut = central::rotatingSurfaceUt(spin, b.state.r, b.state.theta, omega);
    const kerr::Vec4 us = {ut, 0.0, 0.0, omega * ut};
    return std::max(1.0, -kerr::dot(kerr::metricTensor(spin, b.state.r, b.state.theta), u, us));
}

}  // namespace

glm::vec3 displayColor(const Body& body) {
    switch (body.kind) {
        case BodyKind::RockyPlanet: return {0.3f, 0.55f, 0.9f};
        case BodyKind::GasGiant: return {0.85f, 0.7f, 0.5f};
        case BodyKind::MiniBlackHole: return {0.6f, 0.45f, 1.0f};
        case BodyKind::Moon:
            switch (body.palette) {
                case Palette::Io: return {0.9f, 0.8f, 0.35f};
                case Palette::Europa: return {0.85f, 0.8f, 0.7f};
                case Palette::Venus: return {0.95f, 0.88f, 0.65f};
                case Palette::Satellite: return {0.85f, 0.85f, 0.9f};
                default: return {0.62f, 0.6f, 0.57f};
            }
        default: return blackbody::hue(body.temperature);
    }
}

const char* kindName(BodyKind kind) {
    switch (kind) {
        case BodyKind::Star: return "Star";
        case BodyKind::RockyPlanet: return "Planet";
        case BodyKind::GasGiant: return "Gas giant";
        case BodyKind::Pulsar: return "Pulsar";
        case BodyKind::TidalStar: return "Star (tidal)";
        case BodyKind::WhiteDwarf: return "White dwarf";
        case BodyKind::NeutronStar: return "Neutron star";
        case BodyKind::MiniBlackHole: return "Small black hole";
        case BodyKind::Moon: return "Moon, rock or satellite";
    }
    return "?";
}

double pointSourceRadius(double radius, double distance, double pixelAngle, double& dim) {
    const double least = kPointRadiusPixels * pixelAngle * distance;
    dim = 1.0;
    if (radius >= least || least <= 0.0) return radius;
    dim = (radius / least) * (radius / least);
    return least;
}

void History::record(int slot, long long k, const Sample& sample) {
    const int col = column(k);
    const size_t row = static_cast<size_t>(kRows * slot);
    image[row * kColumns + col] = sample.position;
    image[(row + 1) * kColumns + col] = sample.u;
    image[(row + 2) * kColumns + col] = sample.shapeA;
    image[(row + 3) * kColumns + col] = sample.shapeB;
}

namespace {
using central::Host;
using central::hostBit;
constexpr unsigned kAnyHost = 0x1FFu;
constexpr unsigned kHole = hostBit(Host::BlackHole);
}  // namespace

const std::vector<PresetInfo>& presetList() {
    // The orbit demonstrations work around anything (outside the surface, gravity is the
    // same for any body of the same mass); around a star or planet they carry a body
    // that could really be there. The black hole's objects (a star a few kilometres from
    // a neutron star would long since have been torn to gas) are offered around black
    // holes only; stars and planets get what really orbits them.
    static const std::vector<PresetInfo> list = {
        {Preset::IscoOrbit, "ISCO orbit",
         "A circular orbit at the innermost stable circular orbit, the closest a body can circle "
         "without falling in. Closer to the hole for faster spins. Around a star larger than its ISCO, "
         "the innermost orbit skims the surface instead.", false, kAnyHost},
        {Preset::StableOrbit, "Circular orbit", "A comfortable circular orbit at r = 10 scene units.", false,
         kAnyHost},
        {Preset::RetrogradeOrbit, "Retrograde",
         "Circles against the center's rotation. Around a spinning black hole, frame dragging pushes the "
         "retrograde ISCO outward.", false, kAnyHost},
        {Preset::Precessing, "Precessing",
         "An eccentric orbit whose periapsis advances every lap, tracing a rosette (like Mercury, "
         "but far stronger near a black hole or neutron star; imperceptible near a star or planet).", false,
         kAnyHost},
        {Preset::ZoomWhirl, "Zoom-whirl",
         "Tuned just above the separatrix: the body zooms out, then whirls several times around "
         "the hole near an unstable circular orbit before zooming out again.", false, kAnyHost},
        {Preset::PolarOrbit, "Polar orbit",
         "Passes over the poles. With spin, frame dragging (the Lense-Thirring effect) twists the "
         "orbital plane around the spin axis.", false, kAnyHost},
        {Preset::Plunge, "Plunge",
         "Too little angular momentum to stay in orbit: the body spirals in (into the horizon, or onto "
         "the surface).", false, kAnyHost},
        {Preset::Flyby, "Flyby", "An unbound body swings past close to the center and is flung away.", false,
         kAnyHost},
        {Preset::Cluster, "Star cluster", "Ten stars on randomly tilted, near-circular orbits.", false, kHole},
        {Preset::RockyPlanet, "Rocky planet",
         "A rotating Earth-like planet, lit by the accretion disk. With the disk off it goes dark.", true, kHole},
        {Preset::GasGiant, "Gas giant", "A banded giant planet on a wide orbit, lit by the disk.", true, kHole},
        {Preset::Pulsar, "Pulsar",
         "A spinning neutron star sweeping two beams around. You see a flash each time a beam "
         "points at you, including in its lensed images, which flash at different times.", true, kHole},
        {Preset::TidalDisruption, "Tidal disruption",
         "A star falls toward the hole. Its shape follows the tides (the affine star model): an egg "
         "on the way in, squeezed flat at closest approach, then stretched without limit into a long "
         "glowing stream (spaghettification). The debris stays bound and winds around the hole. With "
         "physical units on, the star is a Sun at its true size: a hole above the Hills mass (~1e8 "
         "Suns) swallows it whole.", true, kHole},
        {Preset::TidalDisruptionFlare, "Tidal disruption flare",
         "The full event, as in ESO's animation of AT2019qiz: a hot blue star falls in from far away on a "
         "parabolic orbit and is spaghettified. Half its debris escapes in a long tail; the rest swings out "
         "and comes back one orbit later, crosses itself near the hole (relativistic precession) and shocks. "
         "The hole flares up as a small disk forms, and a billowing outflow erupts above and below. Time "
         "scale ~40.", true, kHole},
        {Preset::GrazingStar, "Grazing star",
         "A star that passes just outside its tidal radius: the tides squeeze and stretch it, and it "
         "leaves wobbling, but it survives.", true, kHole},
        {Preset::BinaryStar, "Binary star",
         "Two stars orbiting each other fall toward the hole. The tide pulls the pair apart (the "
         "Hills mechanism): one star is captured on a tight orbit, the other is flung out fast.", true, kHole},
        {Preset::PlanetWithMoon, "Planet with moon",
         "A giant planet with a moon circling it, both orbiting the hole. The moon stays inside the "
         "planet's Hill sphere, where the planet's pull beats the hole's tide.", true, kHole},
        {Preset::GasCloud, "Gas cloud",
         "A diffuse cloud on an eccentric orbit, like the G2 cloud near our galaxy's central "
         "black hole. It is sheared into a noodle as it swings past.", true,
         kHole | hostBit(Host::NeutronStar) | hostBit(Host::WhiteDwarf)},
        {Preset::DebrisRing, "Tilted ring",
         "A ring of gas tilted 30 degrees to the spin. Frame dragging precesses inner orbits "
         "fastest, so the ring slowly twists and warps.", true,
         kHole | hostBit(Host::NeutronStar) | hostBit(Host::WhiteDwarf)},
        {Preset::MiniBlackHole, "Small black hole",
         "A black hole of 2% of the central mass on a tilted orbit. It lenses whatever lies behind it "
         "and swallows anything that strays within its horizon. (No exact spacetime holds two black "
         "holes, so its lensing is the thin-lens approximation: a bend of 4m/b.)", true, kHole},
        {Preset::WhiteDwarfCompanion, "White dwarf",
         "A hot white dwarf, a stellar core the size of the Earth, on an eccentric orbit.", true, kHole},
        {Preset::Moon, "Moon",
         "The Moon at its real size (1,737 km, a quarter of the Earth's) and distance (384,400 km, 60 Earth "
         "radii), on its tilted, slightly eccentric orbit, always showing the same face. It goes round "
         "once in 27.3 days; lit by the Sun, it passes through the Earth's shadow at a lunar eclipse.", true,
         hostBit(Host::Earthlike)},
        {Preset::SpaceStation, "Space station",
         "The International Space Station, 420 km up on an orbit tilted 51.6 degrees: round the Earth in "
         "92 minutes at 7.7 km/s. At 109 m across it is a point of light unless you are right beside it.", true,
         hostBit(Host::Earthlike)},
        {Preset::GeostationarySatellite, "Geostationary satellite",
         "36,000 km up over the equator, where an orbit takes exactly one sidereal day: it hangs over the "
         "same spot as the Earth turns.", true, hostBit(Host::Earthlike)},
        {Preset::AsteroidImpact, "Asteroid impact",
         "A 10 km asteroid like the one that ended the dinosaurs, arriving at 17 km/s and striking at "
         "20 km/s. It leaves a glowing scar.", true, hostBit(Host::Earthlike)},
        {Preset::GalileanMoons, "Galilean moons",
         "Io, Europa, Ganymede and Callisto at their real sizes (Ganymede is larger than Mercury) and "
         "distances (6 to 27 Jupiter radii), lit by the Sun and eclipsed by Jupiter's shadow. With gravity "
         "between bodies on, Io, Europa and Ganymede keep their 1:2:4 resonance.", true,
         hostBit(Host::GasGiant)},
        {Preset::CometImpact, "Comet impact",
         "A 1 km fragment of a comet like Shoemaker-Levy 9, which struck Jupiter in 1994 at 60 km/s.", true,
         hostBit(Host::GasGiant)},
        {Preset::InnerPlanets, "Mercury, Venus, Earth",
         "The inner planets on their real orbits (0.39, 0.72 and 1 AU: 83, 156 and 215 solar radii) at "
         "their real sizes: from beside the Sun they are points of light. Ride along with the Earth to see "
         "the Sun as we do, half a degree across.", true, hostBit(Host::SunLike)},
        {Preset::HotJupiter, "Hot Jupiter",
         "A giant planet like 51 Pegasi b, the first found around a Sun-like star: 0.05 AU out (11 solar "
         "radii), round in 4.2 days, one face always toward the star.", true, hostBit(Host::SunLike)},
        {Preset::SungrazingComet, "Sungrazing comet",
         "A Kreutz sungrazer: a 2 km nucleus falling from far out, skimming 0.2 solar radii above the "
         "surface at 560 km/s through the corona.", true, hostBit(Host::SunLike)},
        {Preset::ProximaPlanets, "Proxima b and d",
         "Proxima Centauri's planets: b, Earth-sized, in the star's habitable zone 0.049 AU out (11 days a "
         "year), and the smaller d at 0.029 AU. Both are thought to keep one face toward the star. Their "
         "surfaces have never been seen.", true, hostBit(Host::RedDwarf)},
        {Preset::GiantInnerPlanets, "Venus and Earth",
         "Venus and the Earth on their orbits when the Sun is a red giant, 1.6 and 2.2 of its radii out: "
         "scorched, under a thousand times today's sunlight. (Mercury is already inside the giant.)", true,
         hostBit(Host::Giant)},
        {Preset::WhiteDwarfPlanet, "Giant planet",
         "WD 1856+534 b: a planet the size of Jupiter, seven times wider than its white dwarf, circling it "
         "every 1.4 days at 3 million km. It survived its star's red-giant phase.", true,
         hostBit(Host::WhiteDwarf)},
        {Preset::WhiteDwarfAsteroid, "Asteroid torn apart",
         "A 50 km asteroid flung toward the white dwarf. Within ~50 white-dwarf radii (its tidal radius, "
         "about the Sun's radius) the tides tear it apart, and its debris winds into a ring: how white "
         "dwarfs get the rocky pollution seen in their spectra.", true, hostBit(Host::WhiteDwarf)},
        {Preset::CataclysmicCompanion, "Red dwarf companion",
         "A cataclysmic variable: a red dwarf (0.15 solar masses, 0.19 solar radii, 15 times wider than the "
         "white dwarf) so close that it fills its Roche lobe, circling in about two hours. Its gas feeds "
         "the disk.", true, hostBit(Host::WhiteDwarf)},
        {Preset::InfallingClump, "Infalling clump",
         "A clump of gas, the debris of an asteroid shredded far out, falls onto the neutron star. It "
         "strikes at half the speed of light and releases a fifth of its rest energy as a flash of X-rays.",
         true, hostBit(Host::NeutronStar)},
    };
    return list;
}

bool presetAvailable(const PresetInfo& preset, central::Host host) {
    return (preset.hosts & central::hostBit(host)) != 0;
}

void BodySystem::clear() {
    bodies.clear();
    departed.clear();
    events.clear();
    pendingDebris.clear();
    tidalEvents.clear();
    impacts.clear();
    accreted = 0.0;
    std::fill(std::begin(history.length), std::end(history.length), 0);
    m_counter = 0;
}

double BodySystem::boundary(double spin) const {
    return central::boundaryRadius(spin, surface);
}

void BodySystem::setScale(double newScale, double spin) {
    scale = newScale;
    retuneHistory(spin);
}

void BodySystem::setFrameTime(double dt, double spin) {
    if (dt <= 0.0) return;
    m_frameTime = m_frameTime > 0.0 ? 0.9 * m_frameTime + 0.1 * dt : dt;
    retuneHistory(spin);
}

void BodySystem::retuneHistory(double spin) {
    // Light takes longer to cross a bigger scene: the history's step follows the scale
    // whenever the two drift apart by more than a factor sqrt(2). Around a star or
    // planet light crosses the scene thousands of times per orbit, so when time runs
    // fast the step also grows to ~1/64 of a frame's worth of motion (re-laid only when
    // off by more than a factor 2), rather than recording millions of samples a frame.
    const double base = 0.5 * scale;
    const double fast = m_frameTime / 64.0;
    const double wanted = std::max(base, fast);
    const double tolerance = fast > base ? 1.0 : 0.5;
    if (std::abs(std::log2(wanted / history.step)) <= tolerance) return;
    history.step = wanted;
    history.newest = static_cast<long long>(std::floor(time / history.step));
    ++history.generation;
    for (const Departed& d : departed) history.length[d.body.slot] = 0;
    departed.clear();
    for (Body& b : bodies) backfill(b, spin);
}

std::string BodySystem::nextName() {
    return "Star " + std::to_string(++m_counter);
}

void BodySystem::logEvent(std::string text) {
    events.push_front(std::move(text));
    while (events.size() > kMaxEvents) events.pop_back();
}

bool BodySystem::push(Body body, double spin) {
    bool used[History::kSlots] = {};
    for (const Body& b : bodies) used[b.slot] = true;
    for (const Departed& d : departed) used[d.body.slot] = true;
    if (std::find(std::begin(used), std::end(used), false) == std::end(used) && !departed.empty()) {
        // Out of slots: let the oldest departed body's light go.
        history.length[departed.front().body.slot] = 0;
        used[departed.front().body.slot] = false;
        departed.erase(departed.begin());
    }
    int slot = -1;
    for (int i = 0; i < History::kSlots && slot < 0; ++i) {
        if (!used[i]) slot = i;
    }
    if (slot < 0) {
        logEvent("Body limit reached (32). Remove some first.");
        return false;
    }
    body.slot = slot;
    if (body.name.empty()) body.name = nextName();
    if (body.seed == 0.0f) body.seed = std::uniform_real_distribution<float>(1.0f, 100.0f)(m_rng);
    body.state.t = time;
    backfill(body, spin);
    recordTrail(body, spin);
    bodies.push_back(std::move(body));
    return true;
}

void BodySystem::backfill(Body& body, double spin) {
    // Run the orbit backwards so the body has a past: with light delay on, its old
    // light is already on its way to the camera and it does not pop into view late.
    // (The past ignores the other bodies' pull and keeps today's shape.)
    Body past = body;
    kerr::State& s = past.state;
    // Followed back to just outside the horizon, or to the surface it would have come from.
    const double rh = surface > 0.0 ? surface : kerr::horizonRadius(spin) + 0.05;
    long long k = static_cast<long long>(std::floor(s.t / history.step));
    History::Sample later = sampleOf(spin, past);
    // Far out (a planet hundreds of stellar radii away) the light now arriving left it
    // up to ~2r ago: lay down at least that much.
    const int wanted = std::clamp(static_cast<int>(std::ceil(2.0 * body.state.r / history.step)), kBackfillSamples,
                                  History::kColumns - 1);
    int n = 0;
    for (int guard = 0; n < wanted && guard < 300000; ++guard) {
        const double tLater = s.t;
        kerr::stepCoordinateTime(spin, past.constants, s, -kerr::suggestedTimeStep(spin, past.constants, s, 0.01));
        kerr::normalizeAngles(s);
        if (!std::isfinite(s.r) || s.r < rh || s.r > 3000.0 * scale) break;
        const History::Sample earlier = sampleOf(spin, past);
        while (n < wanted && k * history.step >= s.t) {
            const float f = static_cast<float>((k * history.step - s.t) / (tLater - s.t));
            history.record(body.slot, k, mix(earlier, later, f));
            ++n;
            --k;
        }
        later = earlier;
    }
    history.length[body.slot] = n;
    history.dirtyRows.push_back(body.slot);
}

void BodySystem::remove(size_t index) {
    if (index >= bodies.size()) return;
    history.length[bodies[index].slot] = 0;
    bodies.erase(bodies.begin() + static_cast<std::ptrdiff_t>(index));
}

void BodySystem::depart(size_t index) {
    if (index >= bodies.size()) return;
    const double lost = bodies[index].state.t;
    departed.push_back({std::move(bodies[index]), lost});
    bodies.erase(bodies.begin() + static_cast<std::ptrdiff_t>(index));
}

bool BodySystem::addFromZamo(double spin, double r, double theta, double phi, const double v[3], BodyKind kind,
                             const std::string& name) {
    if (r < boundary(spin) + 0.2 * scale) {
        logEvent(surface > 0.0 ? "Cannot place a body inside or at the surface."
                               : "Cannot place a body inside or at the horizon.");
        return false;
    }
    if (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] >= 1.0) {
        logEvent("Speed must be below the speed of light.");
        return false;
    }
    Body body;
    body.name = name;
    body.kind = kind;
    kerr::launchFromZamo(spin, r, theta, phi, v, 1.0, body.constants, body.state);
    return push(std::move(body), spin);
}

bool BodySystem::addCircular(double spin, double r, bool prograde, double phase, BodyKind kind,
                             const std::string& name) {
    double s = prograde ? 1.0 : -1.0;
    double sr = std::sqrt(r);
    if (r * sr - 3.0 * sr + 2.0 * s * spin <= 0.0) {
        logEvent(format("No circular orbit exists at r = %.2f (inside the photon orbit).", r));
        return false;
    }
    if (surface > 0.0 && r < surface) {
        logEvent(format("r = %.3g M is inside the surface (R = %.3g M).", r, surface));
        return false;
    }
    Body body;
    body.name = name;
    body.kind = kind;
    body.constants = kerr::circularOrbit(spin, r, prograde);
    body.state = {time, r, kHalfPi, phase, 0.0, 0.0};
    if (r < kerr::iscoRadius(spin, prograde)) {
        logEvent(format("r = %.2f is inside the ISCO (%.2f): the orbit is unstable.", r,
                        kerr::iscoRadius(spin, prograde)));
    }
    return push(std::move(body), spin);
}

bool BodySystem::addEquatorial(double spin, double rp, double ra, bool prograde, double phase, BodyKind kind,
                               const std::string& name) {
    double E, L;
    if (!kerr::equatorialOrbit(spin, rp, ra, prograde, E, L)) {
        logEvent(format("No bound orbit between r = %.2f and %.2f.", rp, ra));
        return false;
    }
    if (surface > 0.0 && rp < surface) logEvent(format("Periapsis %.3g M is below the surface: it will hit it.", rp));
    Body body;
    body.name = name;
    body.kind = kind;
    body.constants = {E, L, 0.0, 1.0};
    body.state = {time, ra, kHalfPi, phase, 0.0, 0.0};  // start at apoapsis
    return push(std::move(body), spin);
}

bool BodySystem::addCompanion(double spin, size_t primary, double mass, double separation, double tiltDeg,
                              BodyKind kind, const std::string& name) {
    // Put a companion beside the primary (outward along r), and give the pair the
    // relative speed of a circular mutual orbit, sqrt((m1 + m2) / d), shared so that
    // their center of mass keeps the primary's old motion (added relativistically).
    // The relative velocity lies along phi, tilted toward theta by tiltDeg.
    Body& p = bodies[primary];
    double vp[3];
    kerr::zamoVelocity(spin, p.constants, p.state, vp);
    const kerr::LocalAxes ax = kerr::localAxes(spin, p.state.r, p.state.theta, p.state.phi);
    const auto x = kerr::toCartesian(spin, p.state.r, p.state.theta, p.state.phi);
    const auto [r, theta, phi] = kerr::fromCartesian(spin, x[0] + separation * ax.r[0], x[1] + separation * ax.r[1],
                                                     x[2] + separation * ax.r[2]);
    const double m1 = p.mass, m2 = mass;
    const double w = std::sqrt((m1 + m2) / separation);
    const double tilt = tiltDeg * kerr::kPi / 180.0;
    const double dir[3] = {0.0, -std::sin(tilt), std::cos(tilt)};
    const double fp = m2 / (m1 + m2), fc = m1 / (m1 + m2);
    const double uc[3] = {fc * w * dir[0], fc * w * dir[1], fc * w * dir[2]};
    const double up[3] = {-fp * w * dir[0], -fp * w * dir[1], -fp * w * dir[2]};
    double vc[3], vpNew[3];
    gravity::addVelocities(vp, uc, vc);
    gravity::addVelocities(vp, up, vpNew);
    const double t = p.state.t;
    kerr::launchFromZamo(spin, p.state.r, p.state.theta, p.state.phi, vpNew, 1.0, p.constants, p.state);
    p.state.t = t;
    // The primary's past (for light delay) was laid down before the change: redo it.
    backfill(p, spin);
    if (!addFromZamo(spin, r, theta, phi, vc, kind, name)) return false;
    bodies.back().mass = mass;
    return true;
}

bool BodySystem::apoapsisVelocity(double spin, double rp, double ra, double tilt, double v[3]) {
    double E, L;
    if (!kerr::equatorialOrbit(spin, rp, ra, true, E, L)) return false;
    const kerr::Metric m = kerr::metricAt(spin, ra, kHalfPi);
    const double gamma = (E - m.omega * L) / std::sqrt(m.sigma * m.delta / m.A);
    const double vt = L / (gamma * std::sqrt(m.A / m.sigma));
    v[0] = 0.0;
    v[1] = -std::sin(tilt) * vt;
    v[2] = std::cos(tilt) * vt;
    return true;
}

bool BodySystem::addRealOrbit(double spin, const std::string& name, BodyKind kind, double radiusKm, double aKm,
                              double e, double inclDeg, double phase, double massRatio, double rotationSeconds) {
    const double a = central::kmToM(aKm, centralMassSolar);
    const double tilt = inclDeg * kerr::kPi / 180.0;
    double v[3], r;
    if (e < 1e-4) {
        r = a;
        const double speed = kerr::circularZamoSpeed(spin, r, true);
        v[0] = 0.0;
        v[1] = -std::sin(tilt) * speed;
        v[2] = std::cos(tilt) * speed;
    } else {
        r = a * (1.0 + e);
        if (!apoapsisVelocity(spin, a * (1.0 - e), r, tilt, v)) {
            logEvent("Could not set up the orbit of " + name + ".");
            return false;
        }
    }
    if (!addFromZamo(spin, r, kHalfPi, phase, v, kind, name)) return false;
    Body& b = bodies.back();
    b.radius = static_cast<float>(central::kmToM(radiusKm, centralMassSolar));
    b.mass = massRatio;
    b.temperature = 0.0f;
    // Locked bodies turn once per orbit: the mean motion a^(-3/2) (Kepler, in M).
    b.spinRate = rotationSeconds == 0.0
                     ? static_cast<float>(1.0 / (a * std::sqrt(a)))
                     : static_cast<float>(2.0 * kerr::kPi * central::timeSeconds(1.0, centralMassSolar) / rotationSeconds);
    return true;
}

void BodySystem::realTestMass(Body& b) const {
    using central::Host;
    const auto km = [&](double radiusKm) { return static_cast<float>(central::kmToM(radiusKm, centralMassSolar)); };
    b.spinRate = 0.0f;
    switch (host) {
        case Host::BlackHole:
            return;
        case Host::NeutronStar:
            // So close to a neutron star only plasma orbits: a blob of the disk's gas (the
            // orbiting "hot spots" of models of kilohertz oscillations), ~1e7 K.
            b.kind = BodyKind::Star;
            b.name = "Plasma blob";
            b.radius = km(0.3);
            b.temperature = static_cast<float>(central::displayTemperature(1e7));
            return;
        case Host::WhiteDwarf:
            // Rock would be vaporized and torn apart this close: a blob of the disk's gas.
            b.kind = BodyKind::Star;
            b.name = "Plasma blob";
            b.radius = km(200.0);
            b.temperature = static_cast<float>(central::displayTemperature(3e4));
            return;
        case Host::Earthlike:
            b.kind = BodyKind::Moon;
            b.palette = Palette::Satellite;
            b.name = "Satellite";
            b.radius = km(0.02);
            return;
        case Host::GasGiant:
            // A small inner moon like Amalthea or Metis.
            b.kind = BodyKind::Moon;
            b.palette = Palette::Moon;
            b.name = "Moonlet";
            b.radius = km(20.0);
            return;
        default:
            // Around a star: a comet's nucleus.
            b.kind = BodyKind::Moon;
            b.palette = Palette::Comet;
            b.name = "Comet";
            b.radius = km(3.0);
            return;
    }
}

void BodySystem::addPreset(double spin, Preset preset) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    // Presets are laid out in scene units: around a big star or planet everything is
    // scaled up with it (radii by S, speeds by 1/sqrt(S), spin rates by S^-1.5), so the
    // same demonstrations work around any object. For black holes S = 1.
    const double S = scale, vS = 1.0 / std::sqrt(scale), tS = scale * std::sqrt(scale);
    const float fS = static_cast<float>(S);
    // Innermost stable circular orbits outside the object: the ISCO, or just above the surface.
    const double clear = surface > 0.0 ? surface + 0.5 * S : 0.0;
    const double iscoPro = std::max(kerr::iscoRadius(spin, true), clear);
    const double iscoRetro = std::max(kerr::iscoRadius(spin, false), clear);
    const size_t before = bodies.size();

    // Applies to the body just added, if any.
    auto style = [&](float temperature, float radius, float spinRate = 0.0f) {
        if (bodies.size() > before) {
            bodies.back().temperature = temperature;
            bodies.back().radius = radius * fS;
            bodies.back().spinRate = static_cast<float>(spinRate / tS);
        }
    };
    // A star whose shape follows the tides; its mass follows from its tidal radius,
    // r_t = R (M / m)^(1/3). Radii in M.
    auto tidalStar = [&](double rp, double ra, double tiltDeg, double radius, double tidalRadius, const char* name) {
        double v[3];
        if (!apoapsisVelocity(spin, rp, ra, tiltDeg * kerr::kPi / 180.0, v)) {
            logEvent("Could not set up the orbit.");
            return;
        }
        if (addFromZamo(spin, ra, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::TidalStar, name)) {
            bodies.back().temperature = 7000.0f;
            bodies.back().radius = static_cast<float>(radius);
            bodies.back().tidalRadius = static_cast<float>(tidalRadius);
            bodies.back().mass = std::pow(radius / tidalRadius, 3.0);
        }
    };

    switch (preset) {
        case Preset::IscoOrbit:
            addCircular(spin, iscoPro * 1.002, true, uni(m_rng) * 6.28, BodyKind::Star, "ISCO orbiter");
            style(9000.0f, 0.4f);
            break;
        case Preset::StableOrbit:
            addCircular(spin, std::max(10.0 * S, iscoPro + 3.0 * S), true, uni(m_rng) * 6.28);
            style(5800.0f, 0.6f);
            break;
        case Preset::RetrogradeOrbit:
            addCircular(spin, std::max(12.0 * S, iscoRetro * 1.3), false, uni(m_rng) * 6.28, BodyKind::Star,
                        "Retrograde star");
            style(4200.0f, 0.6f);
            break;
        case Preset::Precessing:
            addEquatorial(spin, std::max(iscoPro * 1.5, 6.5 * S), 24.0 * S, true, uni(m_rng) * 6.28, BodyKind::Star,
                          "Precessing star");
            style(7000.0f, 0.5f);
            break;
        case Preset::ZoomWhirl: {
            double rp;
            if (kerr::separatrixPeriapsis(spin, 22.0 * S, true, rp) && rp > clear) {
                addEquatorial(spin, rp * (1.0 + 2e-4), 22.0 * S, true, uni(m_rng) * 6.28, BodyKind::Star,
                              "Zoom-whirl star");
                style(12000.0f, 0.4f);
            } else if (surface > 0.0) {
                logEvent("Zoom-whirl orbits whirl near the photon orbit (r ~ 4 M): inside this object.");
            } else {
                logEvent("Could not find a zoom-whirl orbit for this spin.");
            }
            break;
        }
        case Preset::PolarOrbit: {
            double r = 10.0 * S;
            double v[3] = {0.0, -1.0 / std::sqrt(r - 2.0), 0.0};
            addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Star, "Polar orbiter");
            style(5000.0f, 0.55f);
            break;
        }
        case Preset::Plunge: {
            double r = 14.0 * S;
            double v[3] = {0.0, 0.0, 0.45 * kerr::circularZamoSpeed(spin, r, true)};
            addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Star, "Doomed star");
            style(6500.0f, 0.55f);
            break;
        }
        case Preset::Flyby: {
            double v[3] = {-0.27 * vS, 0.0, 0.13 * vS};
            addFromZamo(spin, 40.0 * S, kHalfPi, -0.6, v, BodyKind::Star, "Flyby star");
            style(20000.0f, 0.7f);
            break;
        }
        case Preset::Cluster: {
            for (int i = 0; i < 10; ++i) {
                double r = (7.0 + 15.0 * uni(m_rng)) * S;
                if (r < iscoPro * 1.2) r = iscoPro * 1.2;
                double tilt = (uni(m_rng) * 2.0 - 1.0) * 1.05;  // up to ~60 degrees
                double speed = kerr::circularZamoSpeed(spin, r, true) * (0.95 + 0.1 * uni(m_rng));
                double v[3] = {0.0, -std::sin(tilt) * speed, std::cos(tilt) * speed};
                if (!addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v)) break;
                bodies.back().temperature = static_cast<float>(3000.0 * std::pow(20000.0 / 3000.0, uni(m_rng)));
                bodies.back().radius = static_cast<float>((0.25 + 0.3 * uni(m_rng)) * S);
            }
            break;
        }
        case Preset::RockyPlanet:
            addCircular(spin, std::max(9.0 * S, iscoPro * 1.5), true, uni(m_rng) * 6.28, BodyKind::RockyPlanet,
                        "Planet");
            style(0.0f, 0.75f, 0.06f);
            break;
        case Preset::GasGiant:
            addCircular(spin, 17.0 * S, true, uni(m_rng) * 6.28, BodyKind::GasGiant, "Gas giant");
            style(0.0f, 1.5f, 0.1f);
            break;
        case Preset::Pulsar: {
            double r = 11.0 * S, tilt = 0.25;
            double speed = kerr::circularZamoSpeed(spin, r, true);
            double v[3] = {0.0, -std::sin(tilt) * speed, std::cos(tilt) * speed};
            addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Pulsar, "Pulsar");
            style(30000.0f, 0.3f, 2.0f);  // one turn every ~3 M
            break;
        }
        case Preset::TidalDisruption:
            if (massSolar > 0.0 && surface <= 0.0) {
                // A Sun-like star at its true size around a hole of massSolar Suns:
                // R* = R_sun in M, m* / M = 1 / massSolar, r_t = R* (M / m*)^(1/3).
                const double rStar = central::kSunRadiusKm / (central::kKmPerSolarMass * massSolar);
                const double rt = central::tidalRadius(rStar, 1.0 / massSolar);
                const double rmb = kerr::marginallyBoundRadius(spin, true);
                if (rt < rmb) {
                    // The tides only win inside the hole's reach: it falls in whole.
                    double v[3] = {0.0, 0.0, 0.3 * kerr::circularZamoSpeed(spin, 30.0, true)};
                    if (addFromZamo(spin, 30.0, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Star, "Swallowed star")) {
                        bodies.back().temperature = 5800.0f;
                        bodies.back().radius = static_cast<float>(std::max(rStar, 0.05));
                        bodies.back().mass = 1.0 / massSolar;
                    }
                    logEvent(format("A Sun-like star would be torn apart only inside r_t = %.3g M, within this hole's "
                                    "reach (%.2f M): it is swallowed whole, as by any hole above the Hills mass "
                                    "(1.1e8 Suns).", rt, rmb));
                } else {
                    tidalStar(std::max(rt / 2.3, iscoPro * 1.3), std::max(40.0, 2.5 * rt), 20.0, rStar, rt,
                              "Doomed star");
                    logEvent(format("A Sun-like star at its true size (%.3g M across) is torn apart within r_t = "
                                    "%.3g M.", 2.0 * rStar, rt));
                }
            } else {
                // A bound, eccentric orbit (periapsis 7 M, apoapsis 40 M) tilted 20 degrees,
                // released at apoapsis, plunging to well inside the tidal radius (16 M). Its
                // debris stays bound, winding around the hole on precessing orbits. (A star
                // falling from far away would leave half its debris escaping and the rest
                // returning only after many thousands of M.)
                tidalStar(std::max(7.0 * S, iscoPro * 1.3), 40.0 * S, 20.0, 1.5 * S, 16.0 * S, "Doomed star");
            }
            break;
        case Preset::TidalDisruptionFlare: {
            // A parabolic orbit (E = 1) in the equatorial plane, periapsis 7.5 M: the star
            // (1.5 M across, tidal radius 16 M) plunges to beta = r_t / r_p ~ 2.1. Its debris
            // straddles E = 1 with a spread dE ~ M R* / r_p^2 ~ 0.03: about half unbound,
            // the rest returning after ~400 M (a = M / 2 dE ~ 15 M). R(r_p) = 0 with E = 1
            // and Q = 0 fixes L: (r_p - 2) L^2 + 4 a L - 2 (a^2 + r_p^2) = 0.
            const double rp = std::max(6.0 * S, iscoPro * 1.3), r0 = 70.0 * S;
            const double rStar = 1.5 * S, rt = 16.0 * S;
            const double A = rp - 2.0, B = 4.0 * spin, C = -2.0 * (spin * spin + rp * rp);
            const double L = (-B + std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
            Body body;
            body.name = "Blue star";
            body.kind = BodyKind::TidalStar;
            body.constants = {1.0, L, 0.0, 1.0};
            const double R0 = kerr::radialPotential(spin, body.constants, r0);
            body.state = {time, r0, kHalfPi, uni(m_rng) * 6.28, -std::sqrt(std::max(R0, 0.0)), 0.0};
            body.temperature = 25000.0f;  // a hot, young B star
            body.brightness = 2.5f;
            // Inside, the gas is as hot as the star's own gravity makes it (virial: internal
            // energy (3/7) G m / R* per unit mass for a gamma = 5/3 polytrope), which on the
            // debris's temperature scale (1.3e7 K per c^2, particles_update.comp) is ~3,000 K.
            // Hotter, its pressure would blow the stream apart, as nothing holds it together
            // (the gas has no self-gravity here).
            body.interiorTemperature = 3000.0f;
            body.radius = static_cast<float>(rStar);
            body.tidalRadius = static_cast<float>(rt);
            body.mass = std::pow(rStar / rt, 3.0);
            if (push(std::move(body), spin)) {
                logEvent(format("A blue star falls toward the hole on a parabolic orbit; it will be torn apart inside "
                                "r_t = %.0f M (closest approach %.1f M).", rt, rp));
            }
            break;
        }
        case Preset::GrazingStar:
            // Periapsis 14 M, just outside its tidal radius of 12 M (beta = 0.86): close
            // to the edge of disruption, which comes near beta = 0.9.
            tidalStar(14.0 * S, 45.0 * S, -12.0, 1.3 * S, 12.0 * S, "Grazing star");
            break;
        case Preset::BinaryStar: {
            // The pair's center of mass falls from 55 M to 6.5 M; the pair is 2.5 M
            // across, so the hole's tide beats their mutual pull near periapsis.
            double v[3];
            if (!apoapsisVelocity(spin, std::max(6.5 * S, iscoPro * 1.3), 55.0 * S, 0.0, v) ||
                !addFromZamo(spin, 55.0 * S, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Star, "Binary A")) {
                logEvent("Could not set up the binary's orbit.");
                break;
            }
            const size_t a = bodies.size() - 1;
            bodies[a].temperature = 9500.0f;
            bodies[a].radius = 0.35f * fS;
            bodies[a].mass = 0.02;
            // Separated along r, orbiting each other in a plane across the orbit.
            if (addCompanion(spin, a, 0.02, 2.5 * S, 90.0, BodyKind::Star, "Binary B")) {
                bodies.back().temperature = 4500.0f;
                bodies.back().radius = 0.35f * fS;
            }
            logEvent("A binary star falls toward the hole.");
            break;
        }
        case Preset::PlanetWithMoon: {
            if (!addCircular(spin, 22.0 * S, true, uni(m_rng) * 6.28, BodyKind::GasGiant, "Giant planet")) break;
            const size_t p = bodies.size() - 1;
            bodies[p].radius = 0.6f * fS;
            bodies[p].spinRate = static_cast<float>(0.1 / tS);
            bodies[p].mass = 0.006;
            if (addCompanion(spin, p, 1e-6, 1.1 * S, 0.0, BodyKind::RockyPlanet, "Moon")) {
                bodies.back().radius = 0.2f * fS;
                bodies.back().spinRate = static_cast<float>(0.05 / tS);
            }
            break;
        }
        case Preset::GasCloud: {
            const double r = 35.0 * S;
            double v[3] = {0.0, -0.08 * vS, 0.6 * kerr::circularZamoSpeed(spin, r, true)};
            pendingDebris.push_back({debris::cloud(spin, r, kHalfPi, uni(m_rng) * 6.28, v, 0.8 * S, true, 100000,
                                                   6000.0f, m_rng, surface),
                                     true});
            logEvent("A gas cloud falls toward the hole.");
            break;
        }
        case Preset::DebrisRing:
            pendingDebris.push_back({debris::ring(spin, std::max(9.0 * S, iscoPro * 1.4), 13.0 * S, 30.0, 150000,
                                                  3500.0f, m_rng, surface),
                                     true});
            logEvent("A tilted debris ring forms.");
            break;
        case Preset::MiniBlackHole: {
            // 2% of the central mass, circling at 14 M on an orbit tilted 25 degrees.
            const double r = std::max(14.0 * S, iscoPro * 1.5), tilt = 25.0 * kerr::kPi / 180.0, m = 0.02;
            const double speed = kerr::circularZamoSpeed(spin, r, true);
            double v[3] = {0.0, -std::sin(tilt) * speed, std::cos(tilt) * speed};
            if (addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::MiniBlackHole, "Small black hole")) {
                bodies.back().mass = m;
                bodies.back().temperature = 0.0f;
                bodies.back().radius = static_cast<float>(3.0 * std::sqrt(3.0) * m);  // its shadow
            }
            break;
        }
        case Preset::WhiteDwarfCompanion:
            addEquatorial(spin, std::max(9.0 * S, iscoPro * 1.5), 30.0 * S, true, uni(m_rng) * 6.28,
                          BodyKind::WhiteDwarf, "White dwarf");
            style(25000.0f, 0.25f);
            break;
        default:
            addRealPreset(spin, preset);
            break;
    }
    // Around anything but a black hole, the orbit demonstrations carry a body that could
    // really be there, at its real size.
    if (host != central::Host::BlackHole) {
        for (const PresetInfo& info : presetList()) {
            if (info.preset != preset || info.object) continue;
            for (size_t i = before; i < bodies.size(); ++i) realTestMass(bodies[i]);
        }
    }
}

void BodySystem::addRealPreset(double spin, Preset preset) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double M = centralMassSolar;
    const auto km = [&](double radiusKm) { return central::kmToM(radiusKm, M); };
    const auto phase = [&]() { return uni(m_rng) * 2.0 * kerr::kPi; };
    constexpr double kAu = central::kAuKm, kDay = 86400.0;
    // Masses of the companions, as fractions of the central mass.
    const auto ofSuns = [&](double massSolar) { return massSolar / M; };
    const double earthMass = central::kEarthMassSolar, jupiterMass = 9.548e-4;
    // An airless body with the given surface.
    const auto rocky = [&](const std::string& name, Palette palette, double radiusKm, double aKm, double e,
                           double inclDeg, double massRatio, double rotationSeconds = 0.0) {
        if (addRealOrbit(spin, name, BodyKind::Moon, radiusKm, aKm, e, inclDeg, phase(), massRatio, rotationSeconds)) {
            bodies.back().palette = palette;
        }
    };
    // Something falling in from r0 (km) at speedKmS, mostly straight down (a fraction
    // `sideways` of its speed across), to strike the surface.
    const auto impactor = [&](const std::string& name, Palette palette, double radiusKm, double r0Km,
                              double speedKmS, double sideways) {
        const double v0 = speedKmS / central::kSpeedOfLightKmS;
        const double v[3] = {-std::sqrt(1.0 - sideways * sideways) * v0, 0.0, sideways * v0};
        const double theta = kHalfPi + 0.4 * (uni(m_rng) - 0.5);
        if (addFromZamo(spin, km(r0Km), theta, phase(), v, BodyKind::Moon, name)) {
            bodies.back().palette = palette;
            bodies.back().radius = static_cast<float>(km(radiusKm));
            bodies.back().mass = 1e-12;
        }
    };

    switch (preset) {
        // --- Around the Earth ---------------------------------------------------------
        case Preset::Moon:
            // Its orbit is tilted 5.1 degrees to the ecliptic, so 18-29 degrees to the
            // Earth's equator as it precesses (every 18.6 years); 23.4 here.
            rocky("Moon", Palette::Moon, central::kMoonRadiusKm, 384400.0, 0.0549, 23.4, 0.0123);
            break;
        case Preset::SpaceStation:
            // Kept with one side toward the Earth, like the real station.
            rocky("Space station", Palette::Satellite, 0.055, central::kEarthRadiusKm + 420.0, 0.0, 51.64, 1e-19);
            break;
        case Preset::GeostationarySatellite:
            // r^3 = G M (1 sidereal day / 2 pi)^2: 42,164 km from the center.
            rocky("Geostationary satellite", Palette::Satellite, 0.02, 42164.0, 0.0, 0.0, 1e-21);
            break;
        case Preset::AsteroidImpact:
            // From 6.7 Earth radii at 17.1 km/s it gains the escape energy and strikes at
            // sqrt(17.1^2 + 11.2^2 (1 - 1/6.7)) = 20 km/s. Carbonaceous: dark.
            impactor("Asteroid (10 km)", Palette::Comet, 5.0, 6.7 * central::kEarthRadiusKm, 17.1, 0.12);
            break;

        // --- Around Jupiter ---------------------------------------------------------------
        case Preset::GalileanMoons:
            // Radii, semi-major axes, eccentricities, inclinations to Jupiter's equator
            // (JPL), masses relative to Jupiter. All keep one face toward Jupiter.
            rocky("Io", Palette::Io, 1821.6, 421700.0, 0.0041, 0.05, 4.70e-5);
            rocky("Europa", Palette::Europa, 1560.8, 671034.0, 0.0090, 0.47, 2.53e-5);
            rocky("Ganymede", Palette::Ganymede, 2634.1, 1070412.0, 0.0013, 0.20, 7.80e-5);
            rocky("Callisto", Palette::Callisto, 2410.3, 1882709.0, 0.0074, 0.19, 5.67e-5);
            break;
        case Preset::CometImpact:
            // Falling from 6.7 Jupiter radii at 23.6 km/s, it strikes at ~60 km/s.
            impactor("Comet fragment", Palette::Comet, 0.5, 6.67 * central::kJupiterRadiusKm, 23.6, 0.1);
            break;

        // --- Around a Sun-like star -----------------------------------------------------------
        case Preset::InnerPlanets:
            // Orbits from JPL, inclinations to the Sun's equator; Mercury turns three times
            // every two orbits (58.6 days), Venus backwards every 243 days.
            rocky("Mercury", Palette::Mercury, 2439.7, 0.38710 * kAu, 0.2056, 3.38, ofSuns(1.660e-7), 58.646 * kDay);
            rocky("Venus", Palette::Venus, 6051.8, 0.72333 * kAu, 0.0068, 3.86, ofSuns(2.448e-6), -243.02 * kDay);
            if (addRealOrbit(spin, "Earth", BodyKind::RockyPlanet, central::kEarthRadiusKm, kAu, 0.0167, 7.155,
                             phase(), ofSuns(earthMass), 86164.1)) {
                bodies.back().seed = 3.7f;  // the same continents as the Earth at the center
            }
            break;
        case Preset::HotJupiter:
            // 51 Pegasi b: 0.0527 AU, 4.23 days; its radius is not known (it does not
            // transit), so a typical hot Jupiter's 1.2 Jupiter radii, swollen by the heat.
            addRealOrbit(spin, "Hot Jupiter", BodyKind::GasGiant, 1.2 * central::kJupiterRadiusKm, 0.0527 * kAu, 0.0,
                         0.0, phase(), ofSuns(0.46 * jupiterMass));
            break;
        case Preset::SungrazingComet: {
            // On a near-parabolic orbit (v^2 = 2M/r) with perihelion q = 1.2 R: angular
            // momentum sqrt(2 M q). Kreutz sungrazers are retrograde (inclined 144 degrees).
            const double r0 = km(6.67 * central::kSunRadiusKm), q = km(1.2 * centralRadiusKm);
            const double v2 = 2.0 / r0 * 0.9995, vt = std::sqrt(2.0 * q) / r0;
            const double tilt = 144.0 * kerr::kPi / 180.0;
            const double v[3] = {-std::sqrt(std::max(v2 - vt * vt, 0.0)), -std::sin(tilt) * vt, std::cos(tilt) * vt};
            if (addFromZamo(spin, r0, kHalfPi, phase(), v, BodyKind::Moon, "Sungrazing comet")) {
                bodies.back().palette = Palette::Comet;
                bodies.back().radius = static_cast<float>(km(1.0));
                bodies.back().mass = 1e-19;
            }
            break;
        }

        // --- Around a red dwarf (Proxima Centauri) ---------------------------------------------
        case Preset::ProximaPlanets:
            // b: 0.0486 AU, 11.2 days, at least 1.07 Earth masses; d: 0.0289 AU, 5.1 days,
            // 0.26 Earth masses (Faria et al. 2022). Radii for rocky compositions.
            rocky("Proxima b", Palette::Rocky, 6600.0, 0.04856 * kAu, 0.02, 0.0, ofSuns(1.07 * earthMass));
            rocky("Proxima d", Palette::Rocky, 5100.0, 0.02885 * kAu, 0.04, 0.0, ofSuns(0.26 * earthMass));
            break;

        // --- Around a red giant ------------------------------------------------------------------
        case Preset::GiantInnerPlanets:
            rocky("Venus", Palette::Scorched, 6051.8, 0.72333 * kAu, 0.0068, 3.86, ofSuns(2.448e-6), -243.02 * kDay);
            rocky("Earth", Palette::Scorched, central::kEarthRadiusKm, kAu, 0.0167, 7.155, ofSuns(earthMass), 86164.1);
            break;

        // --- Around a white dwarf ---------------------------------------------------------------------
        case Preset::WhiteDwarfPlanet:
            // WD 1856+534 b: 0.0204 AU, 1.41 days, 0.93 Jupiter radii (Vanderburg et al. 2020).
            addRealOrbit(spin, "WD 1856+534 b", BodyKind::GasGiant, 0.93 * central::kJupiterRadiusKm, 0.0204 * kAu, 0.0,
                         0.0, phase(), ofSuns(jupiterMass));
            break;
        case Preset::CataclysmicCompanion: {
            // A red dwarf of 0.15 solar masses and 0.19 solar radii filling its Roche lobe:
            // separation a = R2 / f(q), Eggleton's f(q) = 0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3))).
            const double m2 = 0.15, r2 = 0.19 * central::kSunRadiusKm;
            const double q = m2 / M, q3 = std::cbrt(q);
            const double f = 0.49 * q3 * q3 / (0.6 * q3 * q3 + std::log(1.0 + q3));
            if (addRealOrbit(spin, "Red dwarf", BodyKind::Star, r2, r2 / f, 0.0, 0.0, phase(), q)) {
                bodies.back().temperature = 3100.0f;
            }
            logEvent("A red dwarf fills its Roche lobe: its gas spills over to feed the white dwarf's disk.");
            break;
        }
        case Preset::WhiteDwarfAsteroid: {
            // A 50 km rubble pile of density 3 g/cm^3: tidal radius r_t = R (M / m)^(1/3) =
            // (3 M / (4 pi rho))^(1/3), about a solar radius. It falls in from just outside
            // it, toward a periapsis of 0.3 r_t, where the tides win.
            const double radiusKm = 50.0;
            const double massGrams = 4.0 / 3.0 * kerr::kPi * std::pow(radiusKm * 1e5, 3.0) * 3.0;
            const double ratio = massGrams / (M * 1.98847e33);
            const double rStar = km(radiusKm), rt = central::tidalRadius(rStar, ratio);
            const double rp = std::max(0.3 * rt, 1.5 * surface), ra = 3.0 * rt, r0 = 1.1 * rt;
            double E, L;
            if (!kerr::equatorialOrbit(spin, rp, ra, true, E, L)) {
                logEvent("Could not set up the orbit.");
                break;
            }
            Body body;
            body.name = "Asteroid";
            body.kind = BodyKind::TidalStar;
            body.constants = {E, L, 0.0, 1.0};
            body.state = {time, r0, kHalfPi, phase(),
                          -std::sqrt(std::max(kerr::radialPotential(spin, body.constants, r0), 0.0)), 0.0};
            // Heated by the white dwarf's light to T sqrt(R / 2d): ~1,500 K at r_t, a dull red.
            body.temperature = static_cast<float>(centralTemperature * std::sqrt(surface / (2.0 * rt)));
            body.radius = static_cast<float>(rStar);
            body.tidalRadius = static_cast<float>(rt);
            body.mass = ratio;
            if (push(std::move(body), spin)) {
                logEvent(format("An asteroid falls toward the white dwarf; it will be torn apart inside r_t = %.3g "
                                "white-dwarf radii. Speed up time to watch.", rt / surface));
            }
            break;
        }

        // --- Around a neutron star -------------------------------------------------------------------
        case Preset::InfallingClump: {
            // Nearly radial infall from 40 M: it strikes at ~0.5 c and releases
            // 1 - sqrt(1 - 2M/R) ~ 19% of its rest energy.
            const double r = 40.0 * scale;
            const double v[3] = {-0.2, 0.0, 0.25 * kerr::circularZamoSpeed(spin, r, true)};
            pendingDebris.push_back({debris::cloud(spin, r, kHalfPi + 0.2, phase(), v, 0.5 * scale, true, 60000,
                                                   8000.0f, m_rng, surface),
                                     true});
            logEvent("A clump of gas, the debris of an asteroid torn apart far out, falls onto the neutron star.");
            break;
        }
        default:
            break;
    }
}

void BodySystem::recordTrail(Body& body, double spin) {
    auto p = kerr::toCartesian(spin, body.state.r, body.state.theta, body.state.phi);
    glm::vec3 pos(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
    if (body.trail.empty() || glm::length(pos - body.trail.back()) > 0.15f * static_cast<float>(scale)) {
        body.trail.push_back(pos);
        if (body.trail.size() > kMaxTrail) body.trail.pop_front();
    }
}

void BodySystem::disrupt(const Body& body, double spin) {
    // The debris keeps the star's shape and its internal motion, dF/dtau X, on top of
    // its orbital velocity. The matrices are on the Kerr-Schild Cartesian axes; the
    // debris is laid out on Boyer-Lindquist Cartesian axes, turned about z by the
    // difference of the two azimuths at the star.
    double v[3];
    kerr::zamoVelocity(spin, body.constants, body.state, v);
    const double turn = -(std::atan2(spin, body.state.r) + kerr::azimuthShift(spin, body.state.r));
    const tidal::Mat3 R = {{{std::cos(turn), -std::sin(turn), 0.0}, {std::sin(turn), std::cos(turn), 0.0},
                            {0.0, 0.0, 1.0}}};
    const tidal::Mat3 F = tidal::multiply(tidal::multiply(R, body.shape.F), tidal::transpose(R));
    const tidal::Mat3 Fdot = tidal::multiply(tidal::multiply(R, body.shape.Fdot), tidal::transpose(R));
    // Expanding adiabatically (gamma = 5/3) has cooled the gas: T ~ (det F)^(-2/3).
    const float temperature =
        (body.interiorTemperature > 0.0f ? body.interiorTemperature : body.temperature) *
        static_cast<float>(std::pow(std::max(tidal::determinant(body.shape.F), 0.05), -2.0 / 3.0));
    std::vector<ParticleInit> gas = debris::deformedStar(spin, body.state.r, body.state.theta, body.state.phi, v,
                                                         body.radius, F, Fdot, 150000, temperature, m_rng, surface);

    // When its debris will come back, and what that will power (around a black hole;
    // a star or planet at the center just catches it).
    if (surface <= 0.0 && !gas.empty()) {
        tde::Event ev = tde::analyze(spin, time, gas);
        const kerr::Derivatives d = kerr::derivatives(spin, body.constants, body.state);
        const auto x = kerr::toCartesian(spin, body.state.r, body.state.theta, body.state.phi);
        const auto u = kerr::cartesianVelocity(spin, body.state.r, body.state.theta, body.state.phi, d.dr, d.dtheta,
                                               d.dphi);
        std::array<double, 3> n = {x[1] * u[2] - x[2] * u[1], x[2] * u[0] - x[0] * u[2], x[0] * u[1] - x[1] * u[0]};
        const double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (len > 0.0) ev.axis = {n[0] / len, n[1] / len, n[2] / len};
        const tde::RadialOrbit orbit = tde::radialOrbit(spin, body.constants, body.state);
        ev.pericenter = orbit.rp > 0.0 ? orbit.rp : body.state.r;
        // Gas keeping the orbit's angular momentum, sqrt(2 M r_p), circles at 2 r_p.
        ev.circularization = std::max(2.0 * ev.pericenter, kerr::iscoRadius(spin, true) * 1.2);
        ev.starMass = body.mass;
        ev.starRadius = body.radius;
        ev.tidalRadius = body.tidalRadius;
        tidalEvents.push_back(std::move(ev));
        const tde::Event& e = tidalEvents.back();
        if (!e.fallback.empty()) {
            logEvent(format("%.0f%% of the debris is bound. It starts coming back in %.0f M", 100.0 * e.fallback.boundFraction,
                            e.fallback.firstReturn) +
                     format("; the fallback peaks in %.0f M, then fades as t^(%.2f).", e.fallback.peakTime,
                            e.fallback.lateSlope));
        }
    }
    pendingDebris.push_back({std::move(gas), false});
    logEvent(body.name + " was torn apart by tidal forces.");
}

void BodySystem::drift(double spin, double target) {
    const double rh = kerr::horizonRadius(spin);
    const double limit = kerr::interiorLimit(spin);
    const bool star = surface > 0.0;
    std::vector<size_t> finished;

    for (size_t i = 0; i < bodies.size(); ++i) {
        Body& b = bodies[i];
        History::Sample before = sampleOf(spin, b);
        bool gone = false;
        for (int guard = 0; b.state.t < target - 1e-12 && guard < 20000; ++guard) {
            const double t0 = b.state.t;
            const double h = std::min(target - t0, kerr::suggestedTimeStep(spin, b.constants, b.state, 0.01));
            const double rate = clockRate(spin, b.constants, b.state);
            kerr::stepCoordinateTime(spin, b.constants, b.state, h);
            kerr::normalizeAngles(b.state);
            if (b.kind == BodyKind::TidalStar && std::isfinite(b.state.r) && b.state.r > (star ? surface : rh)) {
                tidal::advance(b.shape, tidal::tidalTensor(spin, b.constants, b.state),
                               tidal::frequencySquared(b.tidalRadius), h * rate);
                b.maxStretch = std::max(b.maxStretch, static_cast<float>(tidal::stretch(b.shape.F)));
            }

            const History::Sample after = sampleOf(spin, b);
            // Samples older than the ring holds would only be overwritten.
            const long long firstKept = static_cast<long long>(std::floor(target / history.step)) - History::kColumns + 1;
            for (long long k = std::max(static_cast<long long>(std::floor(t0 / history.step)) + 1, firstKept);
                 k * history.step <= b.state.t + 1e-9; ++k) {
                const float f = static_cast<float>((k * history.step - t0) / (b.state.t - t0));
                history.record(b.slot, k, mix(before, after, f));
                history.length[b.slot] = std::min(History::kColumns, history.length[b.slot] + 1);
            }
            before = after;

            if (!std::isfinite(b.state.r)) {
                logEvent(b.name + " was lost (numerical trouble).");
                gone = true;
                break;
            }
            if (star && b.state.r <= surface + 0.5 * b.radius) {
                // It strikes the surface: the energy of the impact, measured by the surface,
                // is (gamma - 1) of its rest mass. The app turns it into a hot spot.
                const double gamma = impactGamma(spin, b, surfaceOmega);
                impacts.push_back({b.state.theta, b.state.phi, b.state.t, gamma - 1.0, b.mass});
                while (impacts.size() > kMaxImpacts) impacts.erase(impacts.begin());
                accreted += b.mass;
                ++impactCount;
                logEvent(b.name + format(" hit the surface at %.3g c, releasing %.2g of its rest energy.",
                                         std::sqrt(std::max(0.0, 1.0 - 1.0 / (gamma * gamma))), gamma - 1.0));
                gone = true;
            } else if (!star && b.state.r < rh && !b.inside) {
                b.inside = true;
                logEvent(b.name + " crossed the event horizon.");
            }
            if (gone) {
                // Already handled.
            } else if (!star && b.state.r < limit) {
                gone = true;  // as deep as anything is followed
            } else if (b.state.r > 3000.0 * scale) {
                logEvent(b.name + " escaped to infinity.");
                gone = true;
            } else if (b.kind == BodyKind::TidalStar && tidal::stretch(b.shape.F) > kDisruptionStretch) {
                disrupt(b, spin);
                gone = true;
            }
            if (gone) break;
            recordTrail(b, spin);
        }
        if (gone) finished.push_back(i);
        else b.state.t = target;
    }
    for (auto it = finished.rbegin(); it != finished.rend(); ++it) depart(*it);
}

void BodySystem::swallowByMiniHoles(double spin) {
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].kind != BodyKind::MiniBlackHole) continue;
        const kerr::State& h = bodies[i].state;
        const auto xh = kerr::toCartesian(spin, h.r, h.theta, h.phi);
        for (size_t j = 0; j < bodies.size(); ++j) {
            if (j == i) continue;
            const Body& other = bodies[j];
            const auto xo = kerr::toCartesian(spin, other.state.r, other.state.theta, other.state.phi);
            const double d = std::sqrt((xo[0] - xh[0]) * (xo[0] - xh[0]) + (xo[1] - xh[1]) * (xo[1] - xh[1]) +
                                       (xo[2] - xh[2]) * (xo[2] - xh[2]));
            // Its horizon (2m, in units of M) plus the other body's size.
            const double reach = 2.0 * bodies[i].mass + 0.5 * other.radius;
            if (d > reach) continue;
            bodies[i].mass += other.mass;
            logEvent(other.name + " was swallowed by " + bodies[i].name + ".");
            depart(j);
            if (j < i) --i;
            return;  // indices have shifted; the rest waits for the next step
        }
    }
}

void BodySystem::update(double spin, double dt, bool gravity) {
    const double target = time + std::max(dt, 0.0);
    std::vector<gravity::Member> members;
    auto gather = [&]() {
        members.clear();
        for (const Body& b : bodies) members.push_back(memberOf(b));
    };
    auto kick = [&](double h) {
        gather();
        gravity::kick(spin, members, h);
        for (size_t i = 0; i < bodies.size(); ++i) {
            bodies[i].constants = members[i].c;
            bodies[i].state = members[i].s;
        }
    };

    gather();
    if (!gravity || !gravity::interacting(members)) {
        drift(spin, target);
    } else {
        // Kick, drift, kick, in steps short enough for the tightest pair.
        const double maxStep = scale * std::sqrt(scale);  // 1 M, in scene time
        double t = time;
        for (int guard = 0; t < target - 1e-12 && guard < 10000; ++guard) {
            gather();
            const double h = std::min(gravity::stepLimit(spin, members, maxStep), target - t);
            kick(0.5 * h);
            drift(spin, t + h);
            kick(0.5 * h);
            t += h;
        }
    }
    for (const Body& b : bodies) {
        if (b.kind == BodyKind::MiniBlackHole) {
            swallowByMiniHoles(spin);
            break;
        }
    }
    time = target;
    history.newest = static_cast<long long>(std::floor(time / history.step));
    // Light from anything in view left it at most this long ago (see fillGpu).
    const double lightWindow = 800.0 * scale;
    for (size_t i = departed.size(); i-- > 0;) {
        if (time - departed[i].lostTime > lightWindow) {
            history.length[departed[i].body.slot] = 0;
            departed.erase(departed.begin() + static_cast<std::ptrdiff_t>(i));
        }
    }
}

void BodySystem::respin(double oldSpin, double newSpin) {
    const double rhNew = kerr::horizonRadius(newSpin), rhOld = kerr::horizonRadius(oldSpin);
    for (size_t i = bodies.size(); i-- > 0;) {
        Body& b = bodies[i];
        const double t = b.state.t;
        if (b.state.r < kerr::interiorLimit(newSpin) || (surface > 0.0 && b.state.r < surface)) {
            remove(i);
            continue;
        }
        if (b.state.r > rhOld + 0.05 && b.state.r > rhNew + 0.05) {
            double v[3];
            kerr::zamoVelocity(oldSpin, b.constants, b.state, v);
            double speed = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if (speed > 0.999) {
                for (double& c : v) c *= 0.999 / speed;
            }
            kerr::launchFromZamo(newSpin, b.state.r, b.state.theta, b.state.phi, v, 1.0, b.constants, b.state);
        } else {
            // At or inside a horizon there is no ZAMO: keep the constants of motion.
            b.state.vr = -std::sqrt(std::max(kerr::radialPotential(newSpin, b.constants, b.state.r), 0.0));
            b.state.vtheta = std::copysign(
                std::sqrt(std::max(kerr::polarPotential(newSpin, b.constants, b.state.theta), 0.0)), b.state.vtheta);
        }
        b.state.t = t;
        b.inside = b.state.r < rhNew;
    }
}

double BodySystem::reach() const {
    double r = 0.0;
    for (const Body& b : bodies) r = std::max(r, b.state.r + b.radius * b.maxStretch);
    return r;
}

void BodySystem::fillGpu(double spin, bool lightDelay, double cameraR, const glm::dvec3& cameraPosition,
                         double pixelAngle, GpuBodies& out) const {
    out = {};
    // Everything goes out in render units: lengths and times in units of `scale` M.
    const double inv = 1.0 / scale;
    const float fInv = static_cast<float>(inv), fScale = static_cast<float>(scale);
    float shellMin = 1e30f, shellMax = 0.0f;
    int n = 0;
    auto add = [&](const Body& b, double lostAgo) {
        if (n >= GpuBodies::kMax) return;
        const History::Sample now = sampleOf(spin, b);
        GpuBody& g = out.bodies[n++];
        const bool hole = b.kind == BodyKind::MiniBlackHole;
        // Too small to resolve (a planet seen from beside its star): a point of the same
        // total light.
        double dim = 1.0;
        const double distance = glm::length(glm::dvec3(now.position) - cameraPosition);
        const double radius = hole ? b.radius : pointSourceRadius(b.radius, distance, pixelAngle, dim);
        if (radius > b.radius) dim = std::max(dim, kPointLeastBrightness);
        g.positionRadius = glm::vec4(glm::vec3(now.position) * fInv, static_cast<float>(radius * inv));
        // u^t and u^r are dimensionless; u^theta and u^phi are rates per unit length.
        g.fourVelocity = glm::vec4(now.u.x, now.u.y, now.u.z * fScale, now.u.w * fScale);
        // A small black hole carries its mass (render units) where stars keep their
        // temperature; an airless body carries its palette in the seed's hundreds.
        const float seed = b.kind == BodyKind::Moon
                               ? 100.0f * static_cast<float>(b.palette) + std::fmod(b.seed, 100.0f)
                               : b.seed;
        g.params = glm::vec4(hole ? static_cast<float>(b.mass * inv) : b.temperature,
                             b.brightness * static_cast<float>(dim), static_cast<float>(b.kind), seed);
        // Light from it now reaching the camera left it at most about this long ago.
        const int window = std::min(History::kColumns - 1,
                                    static_cast<int>(std::ceil(((cameraR + b.state.r) * 1.2 + 200.0 * scale) /
                                                               history.step)));
        const double angle = std::fmod(b.spinRate * b.state.t, 2.0 * kerr::kPi);
        g.motion = glm::vec4(static_cast<float>(angle), b.spinRate * fScale, static_cast<float>(b.slot),
                             static_cast<float>(history.length[b.slot]));
        g.shapeA = now.shapeA;
        g.shapeB = now.shapeB;

        float rMin = static_cast<float>(b.state.r), rMax = rMin;
        float fastest = 0.0f;  // coordinate speed of the center, for the ray tracer's culling
        if (lightDelay) {
            const int len = std::min(history.length[b.slot], window);
            const glm::vec4* row = &history.image[static_cast<size_t>(History::kRows * b.slot) * History::kColumns];
            for (int j = 0; j < len; ++j) {
                const glm::vec4& sample = row[column(history.newest - j)];
                rMin = std::min(rMin, sample.w);
                rMax = std::max(rMax, sample.w);
                if (j > 0) {
                    const glm::vec4& later = row[column(history.newest - j + 1)];
                    fastest = std::max(fastest, glm::length(glm::vec3(later) - glm::vec3(sample)));
                }
            }
            fastest /= static_cast<float>(history.step);
            const float sinceNewest = static_cast<float>(time - history.newest * history.step);
            if (sinceNewest > 1e-6f * fScale && len > 0) {
                fastest = std::max(fastest, glm::length(glm::vec3(now.position) - glm::vec3(row[column(history.newest)])) /
                                                sinceNewest);
            }
        }
        // Light-travel delay changes a little faster than distance along a ray; keep a margin.
        g.shapeB.w = 1.2f * fastest + 0.01f;
        // The Cartesian mapping stretches distances by up to |a| near the hole, and a
        // tidally deformed star reaches out to its longest axis. A small black hole bends
        // rays passing well outside its shadow.
        const double reach = hole ? kLensReach * b.mass : radius * b.maxStretch;
        const double margin = reach + spin + 0.5 * scale;
        g.range = glm::vec4(static_cast<float>((rMin - margin) * inv), static_cast<float>((rMax + margin) * inv),
                            static_cast<float>(lostAgo * inv), static_cast<float>(reach * 1.01 * inv));
        shellMin = std::min(shellMin, g.range.x);
        shellMax = std::max(shellMax, g.range.y);
    };
    for (const Body& b : bodies) add(b, -1.0);
    // Gone bodies are seen only through their light still in flight.
    if (lightDelay) {
        for (const Departed& d : departed) add(d.body, time - d.lostTime);
    }
    out.info = glm::ivec4(n, column(history.newest), lightDelay ? 1 : 0, History::kColumns);
    out.shell = glm::vec4(std::max(shellMin, 0.0f), shellMax,
                          static_cast<float>((time - history.newest * history.step) / history.step),
                          static_cast<float>(history.step * inv));
}
