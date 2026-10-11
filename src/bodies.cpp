#include "bodies.h"

#include "blackbody.h"
#include "gravity.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {

constexpr double kHalfPi = 0.5 * kerr::kPi;
constexpr size_t kMaxTrail = 3000;
constexpr size_t kMaxEvents = 8;
constexpr int kBackfillSamples = 1400;  // 700 M of past orbit for each new body
// A tidally stretched star is treated as torn apart once its longest axis reaches
// this many radii: from there on its pieces are effectively unbound.
constexpr double kDisruptionStretch = 2.5;
// Gravity between bodies is softened over this fraction of their radii.
constexpr double kSoftening = 0.1;

std::string format(const char* fmt, double a, double b = 0.0) {
    char buf[160];
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

}  // namespace

glm::vec3 displayColor(const Body& body) {
    switch (body.kind) {
        case BodyKind::RockyPlanet: return {0.3f, 0.55f, 0.9f};
        case BodyKind::GasGiant: return {0.85f, 0.7f, 0.5f};
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
    }
    return "?";
}

void History::record(int slot, long long k, const Sample& sample) {
    const int col = column(k);
    const size_t row = static_cast<size_t>(kRows * slot);
    image[row * kColumns + col] = sample.position;
    image[(row + 1) * kColumns + col] = sample.u;
    image[(row + 2) * kColumns + col] = sample.shapeA;
    image[(row + 3) * kColumns + col] = sample.shapeB;
}

const std::vector<PresetInfo>& presetList() {
    static const std::vector<PresetInfo> list = {
        {Preset::IscoOrbit, "ISCO orbit",
         "A circular orbit at the innermost stable circular orbit, the closest a body can circle "
         "without falling in. Closer to the hole for faster spins.", false},
        {Preset::StableOrbit, "Circular orbit", "A comfortable circular orbit at r = 10 M.", false},
        {Preset::RetrogradeOrbit, "Retrograde",
         "Circles against the hole's rotation. Frame dragging pushes the retrograde ISCO outward.", false},
        {Preset::Precessing, "Precessing",
         "An eccentric orbit whose periapsis advances every lap, tracing a rosette (like Mercury, "
         "but far stronger).", false},
        {Preset::ZoomWhirl, "Zoom-whirl",
         "Tuned just above the separatrix: the body zooms out, then whirls several times around "
         "the hole near an unstable circular orbit before zooming out again.", false},
        {Preset::PolarOrbit, "Polar orbit",
         "Passes over the poles. With spin, frame dragging (the Lense-Thirring effect) twists the "
         "orbital plane around the spin axis.", false},
        {Preset::Plunge, "Plunge",
         "Too little angular momentum to stay in orbit: the body spirals in and fades as its "
         "light is redshifted toward the horizon.", false},
        {Preset::Flyby, "Flyby", "An unbound star swings past close to the hole and is flung away.", false},
        {Preset::Cluster, "Star cluster", "Ten stars on randomly tilted, near-circular orbits.", false},
        {Preset::RockyPlanet, "Rocky planet",
         "A rotating Earth-like planet, lit by the accretion disk. With the disk off it goes dark.", true},
        {Preset::GasGiant, "Gas giant", "A banded giant planet on a wide orbit, lit by the disk.", true},
        {Preset::Pulsar, "Pulsar",
         "A spinning neutron star sweeping two beams around. You see a flash each time a beam "
         "points at you, including in its lensed images, which flash at different times.", true},
        {Preset::TidalDisruption, "Tidal disruption",
         "A star falls toward the hole. Its shape follows the tides (the affine star model): an egg "
         "on the way in, squeezed flat at closest approach, then stretched without limit into a long "
         "glowing stream (spaghettification). The debris stays bound and winds around the hole.", true},
        {Preset::GrazingStar, "Grazing star",
         "A star that passes just outside its tidal radius: the tides squeeze and stretch it, and it "
         "leaves wobbling, but it survives.", true},
        {Preset::BinaryStar, "Binary star",
         "Two stars orbiting each other fall toward the hole. The tide pulls the pair apart (the "
         "Hills mechanism): one star is captured on a tight orbit, the other is flung out fast.", true},
        {Preset::PlanetWithMoon, "Planet with moon",
         "A giant planet with a moon circling it, both orbiting the hole. The moon stays inside the "
         "planet's Hill sphere, where the planet's pull beats the hole's tide.", true},
        {Preset::GasCloud, "Gas cloud",
         "A diffuse cloud on an eccentric orbit, like the G2 cloud near our galaxy's central "
         "black hole. It is sheared into a noodle as it swings past.", true},
        {Preset::DebrisRing, "Tilted ring",
         "A ring of rubble tilted 30 degrees to the spin. Frame dragging precesses inner orbits "
         "fastest, so the ring slowly twists and warps.", true},
    };
    return list;
}

void BodySystem::clear() {
    bodies.clear();
    departed.clear();
    events.clear();
    pendingDebris.clear();
    std::fill(std::begin(history.length), std::end(history.length), 0);
    m_counter = 0;
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
    const double rh = kerr::horizonRadius(spin);
    long long k = static_cast<long long>(std::floor(s.t / History::kStep));
    History::Sample later = sampleOf(spin, past);
    int n = 0;
    for (int guard = 0; n < kBackfillSamples && guard < 300000; ++guard) {
        const double tLater = s.t;
        kerr::stepCoordinateTime(spin, past.constants, s, -kerr::suggestedTimeStep(spin, past.constants, s, 0.01));
        kerr::normalizeAngles(s);
        if (!std::isfinite(s.r) || s.r < rh + 0.05 || s.r > 3000.0) break;
        const History::Sample earlier = sampleOf(spin, past);
        while (n < kBackfillSamples && k * History::kStep >= s.t) {
            const float f = static_cast<float>((k * History::kStep - s.t) / (tLater - s.t));
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
    if (r < kerr::horizonRadius(spin) + 0.2) {
        logEvent("Cannot place a body inside or at the horizon.");
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

void BodySystem::addPreset(double spin, Preset preset) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double iscoPro = kerr::iscoRadius(spin, true);
    const double iscoRetro = kerr::iscoRadius(spin, false);
    const size_t before = bodies.size();

    // Applies to the body just added, if any.
    auto style = [&](float temperature, float radius, float spinRate = 0.0f) {
        if (bodies.size() > before) {
            bodies.back().temperature = temperature;
            bodies.back().radius = radius;
            bodies.back().spinRate = spinRate;
        }
    };
    // The ZAMO-frame velocity (tilted out of the equator by tilt) of the equatorial
    // orbit with periapsis rp and apoapsis ra, at apoapsis. False if there is none.
    auto apoapsisVelocity = [&](double rp, double ra, double tilt, double v[3]) {
        double E, L;
        if (!kerr::equatorialOrbit(spin, rp, ra, true, E, L)) return false;
        const kerr::Metric m = kerr::metricAt(spin, ra, kHalfPi);
        const double gamma = (E - m.omega * L) / std::sqrt(m.sigma * m.delta / m.A);
        const double vt = L / (gamma * std::sqrt(m.A / m.sigma));
        v[0] = 0.0;
        v[1] = -std::sin(tilt) * vt;
        v[2] = std::cos(tilt) * vt;
        return true;
    };
    // A star whose shape follows the tides; its mass follows from its tidal radius,
    // r_t = R (M / m)^(1/3).
    auto tidalStar = [&](double rp, double ra, double tiltDeg, float radius, float tidalRadius, const char* name) {
        double v[3];
        if (!apoapsisVelocity(rp, ra, tiltDeg * kerr::kPi / 180.0, v)) {
            logEvent("Could not set up the orbit.");
            return;
        }
        if (addFromZamo(spin, ra, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::TidalStar, name)) {
            style(7000.0f, radius);
            bodies.back().tidalRadius = tidalRadius;
            bodies.back().mass = std::pow(radius / tidalRadius, 3.0);
        }
    };

    switch (preset) {
        case Preset::IscoOrbit:
            addCircular(spin, iscoPro * 1.002, true, uni(m_rng) * 6.28, BodyKind::Star, "ISCO orbiter");
            style(9000.0f, 0.4f);
            break;
        case Preset::StableOrbit:
            addCircular(spin, std::max(10.0, iscoPro + 3.0), true, uni(m_rng) * 6.28);
            style(5800.0f, 0.6f);
            break;
        case Preset::RetrogradeOrbit:
            addCircular(spin, std::max(12.0, iscoRetro * 1.3), false, uni(m_rng) * 6.28, BodyKind::Star,
                        "Retrograde star");
            style(4200.0f, 0.6f);
            break;
        case Preset::Precessing:
            addEquatorial(spin, std::max(iscoPro * 1.5, 6.5), 24.0, true, uni(m_rng) * 6.28, BodyKind::Star,
                          "Precessing star");
            style(7000.0f, 0.5f);
            break;
        case Preset::ZoomWhirl: {
            double rp;
            if (kerr::separatrixPeriapsis(spin, 22.0, true, rp)) {
                addEquatorial(spin, rp * (1.0 + 2e-4), 22.0, true, uni(m_rng) * 6.28, BodyKind::Star,
                              "Zoom-whirl star");
                style(12000.0f, 0.4f);
            } else {
                logEvent("Could not find a zoom-whirl orbit for this spin.");
            }
            break;
        }
        case Preset::PolarOrbit: {
            double r = 10.0;
            double v[3] = {0.0, -1.0 / std::sqrt(r - 2.0), 0.0};
            addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Star, "Polar orbiter");
            style(5000.0f, 0.55f);
            break;
        }
        case Preset::Plunge: {
            double r = 14.0;
            double v[3] = {0.0, 0.0, 0.45 * kerr::circularZamoSpeed(spin, r, true)};
            addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Star, "Doomed star");
            style(6500.0f, 0.55f);
            break;
        }
        case Preset::Flyby: {
            double v[3] = {-0.27, 0.0, 0.13};
            addFromZamo(spin, 40.0, kHalfPi, -0.6, v, BodyKind::Star, "Flyby star");
            style(20000.0f, 0.7f);
            break;
        }
        case Preset::Cluster: {
            for (int i = 0; i < 10; ++i) {
                double r = 7.0 + 15.0 * uni(m_rng);
                if (r < iscoPro * 1.2) r = iscoPro * 1.2;
                double tilt = (uni(m_rng) * 2.0 - 1.0) * 1.05;  // up to ~60 degrees
                double speed = kerr::circularZamoSpeed(spin, r, true) * (0.95 + 0.1 * uni(m_rng));
                double v[3] = {0.0, -std::sin(tilt) * speed, std::cos(tilt) * speed};
                if (!addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v)) break;
                bodies.back().temperature = static_cast<float>(3000.0 * std::pow(20000.0 / 3000.0, uni(m_rng)));
                bodies.back().radius = static_cast<float>(0.25 + 0.3 * uni(m_rng));
            }
            break;
        }
        case Preset::RockyPlanet:
            addCircular(spin, std::max(9.0, iscoPro * 1.5), true, uni(m_rng) * 6.28, BodyKind::RockyPlanet,
                        "Planet");
            style(0.0f, 0.75f, 0.06f);
            break;
        case Preset::GasGiant:
            addCircular(spin, 17.0, true, uni(m_rng) * 6.28, BodyKind::GasGiant, "Gas giant");
            style(0.0f, 1.5f, 0.1f);
            break;
        case Preset::Pulsar: {
            double r = 11.0, tilt = 0.25;
            double speed = kerr::circularZamoSpeed(spin, r, true);
            double v[3] = {0.0, -std::sin(tilt) * speed, std::cos(tilt) * speed};
            addFromZamo(spin, r, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Pulsar, "Pulsar");
            style(30000.0f, 0.3f, 2.0f);  // one turn every ~3 M
            break;
        }
        case Preset::TidalDisruption:
            // A bound, eccentric orbit (periapsis 7 M, apoapsis 40 M) tilted 20 degrees,
            // released at apoapsis, plunging to well inside the tidal radius (16 M). Its
            // debris stays bound, winding around the hole on precessing orbits. (A star
            // falling from far away would leave half its debris escaping and the rest
            // returning only after many thousands of M.)
            tidalStar(std::max(7.0, iscoPro * 1.3), 40.0, 20.0, 1.5f, 16.0f, "Doomed star");
            break;
        case Preset::GrazingStar:
            // Periapsis 14 M, just outside its tidal radius of 12 M (beta = 0.86): close
            // to the edge of disruption, which comes near beta = 0.9.
            tidalStar(14.0, 45.0, -12.0, 1.3f, 12.0f, "Grazing star");
            break;
        case Preset::BinaryStar: {
            // The pair's center of mass falls from 55 M to 6.5 M; the pair is 2.5 M
            // across, so the hole's tide beats their mutual pull near periapsis.
            double v[3];
            if (!apoapsisVelocity(std::max(6.5, iscoPro * 1.3), 55.0, 0.0, v) ||
                !addFromZamo(spin, 55.0, kHalfPi, uni(m_rng) * 6.28, v, BodyKind::Star, "Binary A")) {
                logEvent("Could not set up the binary's orbit.");
                break;
            }
            const size_t a = bodies.size() - 1;
            bodies[a].temperature = 9500.0f;
            bodies[a].radius = 0.35f;
            bodies[a].mass = 0.02;
            // Separated along r, orbiting each other in a plane across the orbit.
            if (addCompanion(spin, a, 0.02, 2.5, 90.0, BodyKind::Star, "Binary B")) {
                bodies.back().temperature = 4500.0f;
                bodies.back().radius = 0.35f;
            }
            logEvent("A binary star falls toward the hole.");
            break;
        }
        case Preset::PlanetWithMoon: {
            if (!addCircular(spin, 22.0, true, uni(m_rng) * 6.28, BodyKind::GasGiant, "Giant planet")) break;
            const size_t p = bodies.size() - 1;
            bodies[p].radius = 0.6f;
            bodies[p].spinRate = 0.1f;
            bodies[p].mass = 0.006;
            if (addCompanion(spin, p, 1e-6, 1.1, 0.0, BodyKind::RockyPlanet, "Moon")) {
                bodies.back().radius = 0.2f;
                bodies.back().spinRate = 0.05f;
            }
            break;
        }
        case Preset::GasCloud: {
            const double r = 35.0;
            double v[3] = {0.0, -0.08, 0.6 * kerr::circularZamoSpeed(spin, r, true)};
            pendingDebris.push_back(
                {debris::cloud(spin, r, kHalfPi, uni(m_rng) * 6.28, v, 0.8, true, 100000, 6000.0f, m_rng), true});
            logEvent("A gas cloud falls toward the hole.");
            break;
        }
        case Preset::DebrisRing:
            pendingDebris.push_back({debris::ring(spin, std::max(9.0, iscoPro * 1.4), 13.0, 30.0, 150000, 3500.0f, m_rng), true});
            logEvent("A tilted debris ring forms.");
            break;
    }
}

void BodySystem::recordTrail(Body& body, double spin) {
    auto p = kerr::toCartesian(spin, body.state.r, body.state.theta, body.state.phi);
    glm::vec3 pos(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
    if (body.trail.empty() || glm::length(pos - body.trail.back()) > 0.15f) {
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
        body.temperature * static_cast<float>(std::pow(std::max(tidal::determinant(body.shape.F), 0.05), -2.0 / 3.0));
    pendingDebris.push_back({debris::deformedStar(spin, body.state.r, body.state.theta, body.state.phi, v, body.radius,
                                                  F, Fdot, 150000, temperature, m_rng),
                             false});
    logEvent(body.name + " was torn apart by tidal forces.");
}

void BodySystem::drift(double spin, double target) {
    const double rh = kerr::horizonRadius(spin);
    const double limit = kerr::interiorLimit(spin);
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
            if (b.kind == BodyKind::TidalStar && std::isfinite(b.state.r) && b.state.r > rh) {
                tidal::advance(b.shape, tidal::tidalTensor(spin, b.constants, b.state),
                               tidal::frequencySquared(b.tidalRadius), h * rate);
                b.maxStretch = std::max(b.maxStretch, static_cast<float>(tidal::stretch(b.shape.F)));
            }

            const History::Sample after = sampleOf(spin, b);
            for (long long k = static_cast<long long>(std::floor(t0 / History::kStep)) + 1;
                 k * History::kStep <= b.state.t + 1e-9; ++k) {
                const float f = static_cast<float>((k * History::kStep - t0) / (b.state.t - t0));
                history.record(b.slot, k, mix(before, after, f));
                history.length[b.slot] = std::min(History::kColumns, history.length[b.slot] + 1);
            }
            before = after;

            if (!std::isfinite(b.state.r)) {
                logEvent(b.name + " was lost (numerical trouble).");
                gone = true;
                break;
            }
            if (b.state.r < rh && !b.inside) {
                b.inside = true;
                logEvent(b.name + " crossed the event horizon.");
            }
            if (b.state.r < limit) {
                gone = true;  // as deep as anything is followed
            } else if (b.state.r > 3000.0) {
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
        double t = time;
        for (int guard = 0; t < target - 1e-12 && guard < 10000; ++guard) {
            gather();
            const double h = std::min(gravity::stepLimit(spin, members, 1.0), target - t);
            kick(0.5 * h);
            drift(spin, t + h);
            kick(0.5 * h);
            t += h;
        }
    }
    time = target;
    history.newest = static_cast<long long>(std::floor(time / History::kStep));
    // Light from anything in view left it at most this long ago (see fillGpu).
    constexpr double kLightWindow = 800.0;
    for (size_t i = departed.size(); i-- > 0;) {
        if (time - departed[i].lostTime > kLightWindow) {
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
        if (b.state.r < kerr::interiorLimit(newSpin)) {
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

void BodySystem::fillGpu(double spin, bool lightDelay, double cameraR, GpuBodies& out) const {
    out = {};
    // Light from anything we can see left it at most about this long ago.
    const int window = std::min(History::kColumns - 1,
                                static_cast<int>(std::ceil((cameraR * 1.5 + 200.0) / History::kStep)));
    float shellMin = 1e30f, shellMax = 0.0f;
    int n = 0;
    auto add = [&](const Body& b, double lostAgo) {
        if (n >= GpuBodies::kMax) return;
        const History::Sample now = sampleOf(spin, b);
        GpuBody& g = out.bodies[n++];
        g.positionRadius = glm::vec4(glm::vec3(now.position), b.radius);
        g.fourVelocity = now.u;
        g.params = glm::vec4(b.temperature, b.brightness, static_cast<float>(b.kind), b.seed);
        const double angle = std::fmod(b.spinRate * b.state.t, 2.0 * kerr::kPi);
        g.motion = glm::vec4(static_cast<float>(angle), b.spinRate, static_cast<float>(b.slot),
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
            fastest /= static_cast<float>(History::kStep);
            const float sinceNewest = static_cast<float>(time - history.newest * History::kStep);
            if (sinceNewest > 1e-6f && len > 0) {
                fastest = std::max(fastest, glm::length(glm::vec3(now.position) - glm::vec3(row[column(history.newest)])) /
                                                sinceNewest);
            }
        }
        // Light-travel delay changes a little faster than distance along a ray; keep a margin.
        g.shapeB.w = 1.2f * fastest + 0.01f;
        // The Cartesian mapping stretches distances by up to |a| near the hole, and a
        // tidally deformed star reaches out to its longest axis.
        const float margin = b.radius * b.maxStretch + static_cast<float>(spin) + 0.5f;
        g.range = glm::vec4(rMin - margin, rMax + margin, static_cast<float>(lostAgo), b.radius * b.maxStretch * 1.01f);
        shellMin = std::min(shellMin, rMin - margin);
        shellMax = std::max(shellMax, rMax + margin);
    };
    for (const Body& b : bodies) add(b, -1.0);
    // Gone bodies are seen only through their light still in flight.
    if (lightDelay) {
        for (const Departed& d : departed) add(d.body, time - d.lostTime);
    }
    out.info = glm::ivec4(n, column(history.newest), lightDelay ? 1 : 0, History::kColumns);
    out.shell = glm::vec4(std::max(shellMin, 0.0f), shellMax,
                          static_cast<float>((time - history.newest * History::kStep) / History::kStep),
                          static_cast<float>(History::kStep));
}
