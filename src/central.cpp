#include "central.h"

#include "kerr.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>

namespace central {

namespace {

constexpr double kPi = kerr::kPi;

// Complete elliptic integral of the first kind K(m), m = k^2, by the AGM.
double ellipticK(double m) {
    double a = 1.0, b = std::sqrt(1.0 - m);
    for (int i = 0; i < 40 && std::abs(a - b) > 1e-16 * a; ++i) {
        const double an = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = an;
    }
    return kPi / (2.0 * a);
}

// Adaptive Simpson quadrature of f over [a, b].
double simpson(const std::function<double(double)>& f, double a, double b, double fa, double fm, double fb,
               double whole, double tol, int depth) {
    const double m = 0.5 * (a + b);
    const double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
    const double flm = f(lm), frm = f(rm);
    const double left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
    const double right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
    // Never ask for more than double precision can give, or the recursion never ends.
    const double floor = 4e-16 * (std::abs(left) + std::abs(right));
    if (depth <= 0 || std::abs(left + right - whole) <= 15.0 * std::max(tol, floor)) {
        return left + right + (left + right - whole) / 15.0;
    }
    return simpson(f, a, m, fa, flm, fm, left, std::max(0.5 * tol, 1e-17), depth - 1) +
           simpson(f, m, b, fm, frm, fb, right, std::max(0.5 * tol, 1e-17), depth - 1);
}

double integrate(const std::function<double(double)>& f, double a, double b, double tol) {
    const double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
    return simpson(f, a, b, fa, fm, fb, (b - a) / 6.0 * (fa + 4.0 * fm + fb), tol, 40);
}

// Boyer-Lindquist metric components at (r, theta): g_tt, g_tphi, g_phiphi.
void blMetric(double a, double r, double theta, double& gtt, double& gtp, double& gpp) {
    const double s2 = std::sin(theta) * std::sin(theta);
    const double sigma = r * r + a * a * (1.0 - s2);
    gtt = -(1.0 - 2.0 * r / sigma);
    gtp = -2.0 * a * r * s2 / sigma;
    gpp = (r * r + a * a + 2.0 * a * a * r * s2 / sigma) * s2;
}

std::string printf1(const char* fmt, double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), fmt, v);
    return buf;
}

Object makeBody(Kind kind, const char* name, double massSolar, double radiusKm, double chi, double inertia,
                double temperature, double field = 0.0, double tiltDeg = 30.0) {
    Object o;
    o.kind = kind;
    o.name = name;
    o.massSolar = massSolar;
    o.radius = radiusKm / (kKmPerSolarMass * massSolar);
    o.chi = chi;
    o.inertia = inertia;
    o.surfaceTemperature = temperature;
    o.magneticField = field;
    o.magneticTiltDeg = tiltDeg;
    return o;
}

Object makeHole(const char* name, double massSolar, double chi) {
    Object o;
    o.kind = Kind::BlackHole;
    o.name = name;
    o.massSolar = massSolar;
    o.chi = chi;
    return o;
}

}  // namespace

// --- Physical scale --------------------------------------------------------------

double lengthKm(double lengthM, double massSolar) {
    return lengthM * kKmPerSolarMass * massSolar;
}

double kmToM(double km, double massSolar) {
    return km / (kKmPerSolarMass * massSolar);
}

double timeSeconds(double timeM, double massSolar) {
    return timeM * kSecondsPerSolarMass * massSolar;
}

double accelerationSI(double accelerationGeometric, double massSolar) {
    const double c = kSpeedOfLightKmS * 1e3;
    return accelerationGeometric * c * c / (kKmPerSolarMass * 1e3 * massSolar);
}

std::string formatLength(double lengthM, double massSolar) {
    const double km = std::abs(lengthKm(lengthM, massSolar));
    if (km < 1.0) return printf1("%.3g m", km * 1e3);
    if (km < 1e3) return printf1("%.3g km", km);
    if (km < 1e7) {
        // Four significant digits with thousands separators: "695,700 km".
        const double unit = std::pow(10.0, std::floor(std::log10(km)) - 3.0);
        const long long n = std::llround(std::round(km / unit) * unit);
        std::string digits = std::to_string(n), grouped;
        for (size_t i = 0; i < digits.size(); ++i) {
            if (i > 0 && (digits.size() - i) % 3 == 0) grouped += ',';
            grouped += digits[i];
        }
        return grouped + " km";
    }
    if (km < 1e4 * kAuKm) return printf1("%.3g AU", km / kAuKm);
    return printf1("%.3g ly", km / kLightYearKm);
}

std::string formatTime(double timeM, double massSolar) {
    const double s = std::abs(timeSeconds(timeM, massSolar));
    if (s < 1e-3) return printf1("%.3g us", s * 1e6);
    if (s < 1.0) return printf1("%.3g ms", s * 1e3);
    if (s < 120.0) return printf1("%.3g s", s);
    if (s < 7200.0) return printf1("%.3g min", s / 60.0);
    if (s < 2.0 * 86400.0) return printf1("%.3g h", s / 3600.0);
    if (s < 2.0 * 3.15576e7) return printf1("%.3g days", s / 86400.0);
    return printf1("%.3g yr", s / 3.15576e7);
}

// --- The object --------------------------------------------------------------------

const char* kindName(Kind kind) {
    switch (kind) {
        case Kind::BlackHole: return "Black hole";
        case Kind::NeutronStar: return "Neutron star";
        case Kind::WhiteDwarf: return "White dwarf";
        case Kind::Star: return "Star";
        case Kind::RockyPlanet: return "Rocky planet";
        case Kind::GasGiant: return "Gas giant";
    }
    return "?";
}

double Object::metricSpin() const {
    if (hasHorizon()) return std::clamp(chi, 0.0, 0.999);
    return std::clamp(chi, 0.0, 0.99);
}

double Object::boundaryRadius() const {
    return hasHorizon() ? kerr::horizonRadius(metricSpin()) : radius;
}

double Object::surfaceOmega() const {
    if (hasHorizon()) {
        const double a = metricSpin();
        return a / (2.0 * kerr::horizonRadius(a));
    }
    // J = I Omega with J = chi M^2 (M = 1), never faster than mass shedding allows.
    const double omega = chi / (inertia * radius * radius);
    return std::min(omega, 0.65 / std::sqrt(radius * radius * radius));
}

double Object::maxChi() const {
    if (hasHorizon()) return 0.999;
    // Omega_max ~ 0.65 sqrt(M / R^3) for rigid relativistic rotation, chi = I Omega / M^2.
    return inertia * 0.65 * std::sqrt(radius);
}

double Object::sceneScale() const {
    return hasHorizon() ? 1.0 : std::max(1.0, radius / 6.0);
}

double boundaryRadius(double a, double surface) {
    return surface > 0.0 ? surface : kerr::horizonRadius(a);
}

Host hostOf(const Object& o) {
    switch (o.kind) {
        case Kind::BlackHole: return Host::BlackHole;
        case Kind::NeutronStar: return Host::NeutronStar;
        case Kind::WhiteDwarf: return Host::WhiteDwarf;
        case Kind::RockyPlanet: return Host::Earthlike;
        case Kind::GasGiant: return Host::GasGiant;
        case Kind::Star: break;
    }
    // Ordinary stars by temperature and surface gravity (the same lines that choose how
    // their surfaces are drawn, src/surface.cpp): brown dwarfs below 2,500 K, giants
    // below log g = 3.5 (cgs), red dwarfs below 3,900 K.
    if (o.surfaceTemperature < 2500.0) return Host::BrownDwarf;
    const double rCm = lengthKm(o.radius, o.massSolar) * 1e5;
    const double g = 6.674e-8 * 1.98847e33 * o.massSolar / (rCm * rCm);
    if (std::log10(g) < 3.5) return Host::Giant;
    if (o.surfaceTemperature < 3900.0) return Host::RedDwarf;
    return Host::SunLike;
}

double chiFromPeriod(double periodSeconds, double massSolar, double radiusKm, double inertia) {
    // In units of M: Omega = 2 pi / P times G M / c^3, R = radiusKm / (G M / c^2).
    const double omega = 2.0 * kPi / periodSeconds * timeSeconds(1.0, massSolar);
    const double R = kmToM(radiusKm, massSolar);
    return inertia * R * R * omega;
}

// --- Catalogue ------------------------------------------------------------------

const std::vector<CatalogEntry>& catalog() {
    // Spins are chi = c J / (G M^2) = I Omega / M^2 from the observed rotation periods:
    // the Earth (one sidereal day) ~740, Jupiter (9.9 h) ~520, the Sun (25.4 days) 0.22,
    // the Crab pulsar (33.7 ms) 0.015. Planets far exceed 1; only black holes cannot.
    static const std::vector<CatalogEntry> list = {
        {"Earth", "Our planet: 3e-6 solar masses, 6,371 km, turning once a sidereal day. Its J/M^2 ~ 740 is "
                  "normal for a planet; only black holes are limited to 1. Add the Moon, a space station or a "
                  "geostationary satellite under Objects, at their real sizes and distances.",
         makeBody(Kind::RockyPlanet, "Earth", kEarthMassSolar, kEarthRadiusKm,
                  chiFromPeriod(86164.1, kEarthMassSolar, kEarthRadiusKm, 0.3307), 0.3307, 288.0),
         32.0, false},
        {"Jupiter", "The largest planet: 0.001 solar masses and 11 times the Earth's size, spinning every 9.9 "
                    "hours. Add its four big moons under Objects.",
         makeBody(Kind::GasGiant, "Jupiter", 9.548e-4, kJupiterRadiusKm,
                  chiFromPeriod(35730.0, 9.548e-4, kJupiterRadiusKm, 0.254), 0.254, 165.0),
         32.0, false},
        {"Brown dwarf", "A failed star, 50 Jupiters but no bigger than Jupiter, too light to fuse hydrogen: it glows "
                        "dull red from its heat of formation (~1,500 K) and spins in a few hours.",
         makeBody(Kind::Star, "Brown dwarf", 0.05, 70000.0, chiFromPeriod(5.0 * 3600.0, 0.05, 70000.0, 0.205), 0.205,
                  1500.0),
         26.0, false},
        {"Red dwarf", "Proxima Centauri: 0.12 solar masses, 15% of the Sun's size, a cool (3,000 K) red star of "
                      "the most common kind, turning once in 83 days. Its planets are under Objects.",
         makeBody(Kind::Star, "Red dwarf", 0.1221, 107280.0, chiFromPeriod(83.0 * 86400.0, 0.1221, 107280.0, 0.205),
                  0.205, 3042.0),
         26.0, false},
        {"Sun", "The Sun: 695,700 km, 5,772 K, turning every 25.4 days (spin 0.22). Its gravitational redshift is "
                "two parts in a million. Mercury, Venus and the Earth are under Objects.",
         makeBody(Kind::Star, "Sun", 1.0, kSunRadiusKm, chiFromPeriod(25.38 * 86400.0, 1.0, kSunRadiusKm, 0.07), 0.07,
                  kSunTemperature),
         26.0, false},
        {"Red giant", "A Sun in old age, swollen to 100 times its size (half an AU, past Mercury's orbit) and cooled "
                      "to 4,000 K. Its spin is the Sun's: J/M^2 does not change as it swells.",
         makeBody(Kind::Star, "Red giant", 1.0, 100.0 * kSunRadiusKm, 0.2, 0.1, 4000.0), 26.0, false},
        {"White dwarf", "A typical white dwarf: 0.6 solar masses in a ball 1.4 times the Earth's size (8,700 km), "
                        "held up by degenerate electrons. Here feeding from a disk that reaches its surface, as in a "
                        "dwarf nova, with the disk's real temperature (ultraviolet-hot inside). Raise the magnetic "
                        "field to make an intermediate polar.",
         // The companion's tide truncates the disk at 0.6 a / (1 + q) = 27 radii (Paczynski
         // 1977); the companion fills its Roche lobe 57 radii away. Framed to show both.
         makeBody(Kind::WhiteDwarf, "White dwarf", 0.6, 8700.0, 0.004, 0.205, 15000.0, 0.0), 620.0, true, 27.0, 72.0},
        {"Sirius B", "The companion of Sirius: 1.02 solar masses in 5,840 km (smaller than the Earth), 25,000 K. Its "
                     "gravitational redshift (80 km/s) was measured with the Hubble Space Telescope.",
         makeBody(Kind::WhiteDwarf, "Sirius B", 1.02, 5840.0, 0.002, 0.15, 25000.0), 32.0, false},
        {"Neutron star", "A young neutron star like the Vela pulsar: 1.4 solar masses in 12 km (R = 5.8 M), the size "
                         "of a city, turning 11 times a second with a 3e12 G field tilted 70 degrees from its spin. "
                         "Its two beams sweep round like a lighthouse's (time slowed 28 times), flickering from pulse "
                         "to pulse. Now and then its crust slips (a glitch) and rings at ~22 Hz, shaking the beams.",
         [] {
             // Vela: P = 89.3 ms, B = 3.4e12 G, ~0.7-1e6 K at the surface (Pavlov et al. 2001).
             Object o = makeBody(Kind::NeutronStar, "Neutron star", 1.4, 12.0, chiFromPeriod(0.0893, 1.4, 12.0, 0.35),
                                 0.35, 1.0e6, 3.4e12, 70.0);
             o.pulsar = true;
             o.beacon = true;
             return o;
         }(),
         // 12,950 M per turn: about 2.6 seconds at 5,000 M a second. Seen 55 degrees from
         // the spin axis, nearer to it than the beams (70 degrees): from there each beam
         // goes all the way round on screen. (From outside their cone, the camera would
         // see them swing to and fro like wipers, though the star turns just the same.)
         75.0, false, 0.0, 55.0, 5000.0},
        {"Pulsar", "The Crab pulsar: a neutron star born in the supernova of 1054, turning 30 times a second, with a "
                   "4e12 G field tilted 60 degrees from its spin. Alone: no disk. Its hot polar caps flash each "
                   "time a beam sweeps past you.",
         [] {
             Object o = makeBody(Kind::NeutronStar, "Crab pulsar", 1.4, 12.0, chiFromPeriod(0.0337, 1.4, 12.0, 0.35),
                                 0.35, 1.5e6, 3.8e12, 60.0);
             o.pulsar = true;
             return o;
         }(),
         // 4,900 M per turn: about five seconds.
         70.0, false, 0.0, 0.0, 1000.0},
        {"Millisecond pulsar", "An accreting millisecond X-ray pulsar like SAX J1808.4-3658: spun up by its disk to "
                               "401 turns a second. Its 2e8 G field holds the disk off at ~26 km; inside, the gas "
                               "slides along the field lines onto the magnetic poles, which glow at ~1e7 K.",
         [] {
             Object o = makeBody(Kind::NeutronStar, "Millisecond pulsar", 1.4, 12.0,
                                 chiFromPeriod(1.0 / 401.0, 1.4, 12.0, 0.35), 0.35, 1.0e6, 2e8, 20.0);
             o.pulsar = true;
             return o;
         }(),
         150.0, true, 13.0, 62.0, 30.0, true},
        {"Magnetar", "A neutron star with a 1e15 G field (a thousand times a pulsar's), turning once every 7.5 s. "
                     "Its decaying field keeps it hot (~4e6 K) with a hotter spot; now and then its crust cracks "
                     "and it bursts.",
         [] {
             Object o = makeBody(Kind::NeutronStar, "Magnetar", 1.4, 12.0, chiFromPeriod(7.5, 1.4, 12.0, 0.35), 0.35,
                                 4.0e6, 1e15, 45.0);
             o.magnetar = true;
             return o;
         }(),
         // 1.1e6 M per turn: 11 seconds, 0.7 of the real pace.
         70.0, false, 0.0, 0.0, 1e5},
        {"Stellar black hole", "A 10 solar-mass black hole, like those found in X-ray binaries.",
         makeHole("Stellar black hole", 10.0, 0.7), 32.0, true},
        {"Sgr A*", "The black hole at the center of the Milky Way: 4.3 million Suns. 1 M of time is 21 s.",
         makeHole("Sgr A*", 4.3e6, 0.9), 32.0, true},
        {"M87*", "The first black hole imaged: 6.5 billion Suns. 1 M is 64 AU across; 1 M of time is 9 hours.",
         makeHole("M87*", 6.5e9, 0.9), 32.0, true},
    };
    return list;
}

// --- Sizes ---------------------------------------------------------------------------

std::string sizeComparison(double radiusKm) {
    struct Reference {
        const char* name;
        double diameterKm;
    };
    static const Reference refs[] = {
        {"a city (Manhattan is 22 km long)", 22.0},
        {"the Moon", 2.0 * kMoonRadiusKm},
        {"the Earth", 2.0 * kEarthRadiusKm},
        {"Jupiter", 2.0 * kJupiterRadiusKm},
        {"the Sun", 2.0 * kSunRadiusKm},
        {"Mercury's orbit", 2.0 * 0.3871 * kAuKm},
        {"the Earth's orbit", 2.0 * kAuKm},
    };
    const double d = 2.0 * radiusKm;
    // The reference nearest in size (on a log scale).
    const Reference* best = &refs[0];
    for (const Reference& ref : refs) {
        if (std::abs(std::log(d / ref.diameterKm)) < std::abs(std::log(d / best->diameterKm))) best = &ref;
    }
    std::string out = formatLength(kmToM(d, 1.0), 1.0) + " across: ";
    const double ratio = d / best->diameterKm;
    if (ratio > 0.87 && ratio < 1.15) return out + "about the size of " + best->name;
    return out + printf1("%.2g times the size of ", ratio) + best->name;
}

// --- Light -----------------------------------------------------------------------------

double fluxRelativeToEarth(double radiusKm, double temperature, double distanceKm) {
    const double r = radiusKm / kSunRadiusKm, t = temperature / kSunTemperature, d = kAuKm / distanceKm;
    return r * r * t * t * t * t * d * d;
}

// --- Surface physics --------------------------------------------------------------

double surfaceRedshift(double R) {
    return 1.0 / std::sqrt(1.0 - 2.0 / R) - 1.0;
}

double surfaceGravity(double R) {
    return 1.0 / (R * R * std::sqrt(1.0 - 2.0 / R));
}

double escapeSpeed(double R) {
    return std::sqrt(2.0 / R);
}

double surfaceClockRate(double R) {
    return std::sqrt(1.0 - 2.0 / R);
}

bool photonSphereVisible(double R) {
    return R < 3.0;
}

bool iscoOutside(double a, double R) {
    return R < kerr::iscoRadius(a, true);
}

double rotatingSurfaceUt(double a, double R, double theta, double omega) {
    double gtt, gtp, gpp;
    blMetric(a, R, theta, gtt, gtp, gpp);
    const double norm = -(gtt + 2.0 * omega * gtp + omega * omega * gpp);
    return norm > 0.0 ? 1.0 / std::sqrt(norm) : 0.0;
}

double bendingAngle(double R, double alpha) {
    // psi = integral_0^{1/R} du / sqrt(1/b^2 - u^2 (1 - 2u)), b = R sin(alpha) / sqrt(1 - 2/R).
    // With u = u_R (1 - s^2) the inverse square root at the turning point (alpha = pi/2)
    // becomes finite.
    const double uR = 1.0 / R;
    const double sa = std::sin(alpha);
    if (sa < 1e-12) return 0.0;
    const double fR = uR * uR * (1.0 - 2.0 * uR);  // u^2 (1 - 2u) at the surface
    const double invB2 = fR / (sa * sa);
    if (R <= 3.0 && invB2 <= 1.0 / 27.0) return INFINITY;  // falls back: never escapes
    const double cot2 = alpha >= 0.5 * kPi ? 0.0 : (std::cos(alpha) * std::cos(alpha)) / (sa * sa);
    // 1/b^2 - u^2 (1 - 2u) = f_R cot^2(alpha) + (u_R - u) [u_R + u - 2 (u_R^2 + u_R u + u^2)],
    // factored so that nothing cancels near the turning point; u_R - u = u_R s^2.
    auto f = [&](double s) {
        const double u = uR * (1.0 - s * s);
        const double g = uR + u - 2.0 * (uR * uR + uR * u + u * u);
        if (s <= 0.0) return cot2 > 0.0 ? 0.0 : 2.0 * uR / std::sqrt(uR * g);
        return 2.0 * uR / std::sqrt(fR * cot2 / (s * s) + uR * g);
    };
    return integrate(f, 0.0, 1.0, 1e-13);
}

double bendingAngleBeloborodov(double R, double alpha) {
    const double c = 1.0 - (1.0 - std::cos(alpha)) / (1.0 - 2.0 / R);
    return std::acos(std::clamp(c, -1.0, 1.0));
}

double maxBendingAngle(double R) {
    if (R <= 3.0) return INFINITY;  // light can wind around the star: everything shows
    return bendingAngle(R, 0.5 * kPi);
}

double visibleFraction(double R) {
    const double psi = maxBendingAngle(R);
    return psi >= kPi ? 1.0 : 0.5 * (1.0 - std::cos(psi));
}

double fullVisibilityRadius() {
    double lo = 3.0001, hi = 4.0;  // psi_max(lo) > pi > psi_max(hi)
    for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        (maxBendingAngle(mid) > kPi ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double periapsisAdvance(double p, double e) {
    const double d = p - 6.0 + 2.0 * e;
    return 4.0 * std::sqrt(p / d) * ellipticK(4.0 * e / d) - 2.0 * kPi;
}

// --- Accretion -----------------------------------------------------------------------

double blackHoleEfficiency(double a) {
    return 1.0 - kerr::circularOrbit(a, kerr::iscoRadius(a, true), true).E;
}

SurfaceAccretion surfaceEfficiency(double a, double R, double omega) {
    SurfaceAccretion out;
    out.rInner = std::max(kerr::iscoRadius(a, true), R);
    const double eIn = kerr::circularOrbit(a, out.rInner, true).E;
    // At rest on the equator of the rotating surface: E = -u_t = -u^t (g_tt + Omega g_tphi).
    double gtt, gtp, gpp;
    blMetric(a, R, 0.5 * kPi, gtt, gtp, gpp);
    const double ut = rotatingSurfaceUt(a, R, 0.5 * kPi, omega);
    const double eSurface = -ut * (gtt + omega * gtp);
    out.disk = 1.0 - eIn;
    out.boundaryLayer = eIn - eSurface;
    out.total = 1.0 - eSurface;
    return out;
}

double alfvenRadiusKm(double fieldGauss, double radiusKm, double massSolar, double mdot) {
    constexpr double G = 6.674e-8;          // cgs
    constexpr double sunGrams = 1.98847e33;
    if (fieldGauss <= 0.0 || mdot <= 0.0) return 0.0;
    const double rCm = radiusKm * 1e5;
    const double mu = fieldGauss * rCm * rCm * rCm;
    const double mu2 = mu * mu;
    const double rA = std::pow(mu2 * mu2 / (2.0 * G * massSolar * sunGrams * mdot * mdot), 1.0 / 7.0);
    return rA * 1e-5;
}

double lightCylinderKm(double periodSeconds) {
    return kSpeedOfLightKmS * std::abs(periodSeconds) / (2.0 * kPi);
}

double radioEmissionHeightKm(double periodSeconds) {
    constexpr double heightAtOneSecondKm = 400.0;  // 40 neutron-star radii of 10 km
    const double p = std::abs(periodSeconds);
    return std::min(heightAtOneSecondKm * std::pow(p, 0.3), lightCylinderKm(p));
}

double beamHalfAngle(double periodSeconds) {
    const double sin2 = std::clamp(radioEmissionHeightKm(periodSeconds) / lightCylinderKm(periodSeconds), 0.0, 1.0);
    const double theta = std::asin(std::sqrt(sin2));
    // Past theta = 54.7 deg (tan^2 = 2) the lines' tangents tip past the magnetic equator:
    // the beam then covers the whole hemisphere.
    return std::min(theta + std::atan(0.5 * std::tan(std::min(theta, 0.5 * kPi - 1e-9))), 0.5 * kPi);
}

double beamCoreHalfWidth(double periodSeconds) {
    return 1.22 * kPi / 180.0 / std::sqrt(std::max(periodSeconds, 1e-6));
}

double crustTorsionalFrequency(int l) {
    constexpr double kF2 = 22.0;  // Hz
    return kF2 * std::sqrt((l - 1.0) * (l + 2.0) / 4.0);
}

double rocheLobeFraction(double q) {
    const double q3 = std::cbrt(q);
    return 0.49 * q3 * q3 / (0.6 * q3 * q3 + std::log(1.0 + q3));
}

double innerLagrangeFraction(double q) {
    // The gradient runs from -infinity next to M1 to +infinity next to M2: bisect.
    auto force = [q](double x) { return -1.0 / (x * x) + q / ((1.0 - x) * (1.0 - x)) + (1.0 + q) * x - q; };
    double lo = 1e-6, hi = 1.0 - 1e-6;
    for (int k = 0; k < 200; ++k) {
        const double mid = 0.5 * (lo + hi);
        (force(mid) < 0.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double tidalTruncationFraction(double q) {
    return 0.6 / (1.0 + q);
}

std::vector<std::array<double, 2>> rocheStream(double q, double stop, int count, double& closest) {
    // The restricted three-body problem in the frame turning with the binary (a = 1,
    // G (M1 + M2) = 1, Omega = 1), centered on the center of mass: M1 at -mu2, M2 at mu1.
    const double mu1 = 1.0 / (1.0 + q), mu2 = q / (1.0 + q);
    const double x1 = -mu2, x2 = mu1;
    auto accel = [&](double x, double y, double vx, double vy, double& ax, double& ay) {
        const double dx1 = x - x1, dx2 = x - x2;
        const double r1 = std::sqrt(dx1 * dx1 + y * y), r2 = std::sqrt(dx2 * dx2 + y * y);
        const double k1 = mu1 / (r1 * r1 * r1), k2 = mu2 / (r2 * r2 * r2);
        ax = 2.0 * vy + x - k1 * dx1 - k2 * dx2;
        ay = -2.0 * vx + y - k1 * y - k2 * y;
    };
    // From rest just inside L1 (it is a saddle: the gas falls off it toward M1).
    double x = x1 + innerLagrangeFraction(q) - 1e-4, y = 0.0, vx = 0.0, vy = 0.0;
    std::vector<std::array<double, 2>> path{{x - x1, y}};
    const double dt = 2e-4;
    double r = std::hypot(x - x1, y), last = r;
    closest = r;
    std::vector<std::array<double, 2>> all{{x - x1, y}};
    for (int step = 0; step < 200000; ++step) {
        // Classical RK4 on (x, y, vx, vy).
        double a1x, a1y, a2x, a2y, a3x, a3y, a4x, a4y;
        accel(x, y, vx, vy, a1x, a1y);
        accel(x + 0.5 * dt * vx, y + 0.5 * dt * vy, vx + 0.5 * dt * a1x, vy + 0.5 * dt * a1y, a2x, a2y);
        accel(x + 0.5 * dt * (vx + 0.5 * dt * a1x), y + 0.5 * dt * (vy + 0.5 * dt * a1y), vx + 0.5 * dt * a2x,
              vy + 0.5 * dt * a2y, a3x, a3y);
        accel(x + dt * (vx + 0.5 * dt * a2x), y + dt * (vy + 0.5 * dt * a2y), vx + dt * a3x, vy + dt * a3y, a4x, a4y);
        x += dt * (vx + dt / 6.0 * (a1x + a2x + a3x));
        y += dt * (vy + dt / 6.0 * (a1y + a2y + a3y));
        vx += dt / 6.0 * (a1x + 2.0 * a2x + 2.0 * a3x + a4x);
        vy += dt / 6.0 * (a1y + 2.0 * a2y + 2.0 * a3y + a4y);
        r = std::hypot(x - x1, y);
        all.push_back({x - x1, y});
        closest = std::min(closest, r);
        if (r <= stop || r > last) break;  // reached the disk, or passed its closest approach
        last = r;
    }
    // Resample to `count` points, evenly in time.
    path.clear();
    const int n = std::max(count, 2);
    for (int k = 0; k < n; ++k) {
        const size_t i = static_cast<size_t>(std::lround(static_cast<double>(k) / (n - 1) * (all.size() - 1)));
        path.push_back(all[i]);
    }
    return path;
}

double eddingtonRate(double massSolar) {
    constexpr double c = 2.99792458e10;
    return 1.26e38 * massSolar / (0.1 * c * c);
}

double diskPeakTemperature(double massSolar, double mdot, double rInnerKm) {
    constexpr double G = 6.674e-8, sigma = 5.6704e-5, sunGrams = 1.98847e33;  // cgs
    if (mdot <= 0.0 || rInnerKm <= 0.0) return 0.0;
    const double r = rInnerKm * 1e5;
    const double tStar = std::pow(3.0 * G * massSolar * sunGrams * mdot / (8.0 * kPi * sigma * r * r * r), 0.25);
    // x^(-3/4) (1 - x^(-1/2))^(1/4) at its peak, x = r / r_in = 49/36.
    const double x = 49.0 / 36.0;
    return tStar * std::pow(x, -0.75) * std::pow(1.0 - 1.0 / std::sqrt(x), 0.25);
}

// --- Tidal disruption -----------------------------------------------------------------

double tidalRadius(double rStar, double massRatio) {
    return rStar * std::pow(massRatio, -1.0 / 3.0);
}

double hillsMass(double rStarKm, double mStarSolar) {
    // R* c^2 / (2G) as a mass: G M / c^2 = R* / 2.
    const double m = 0.5 * rStarKm / kKmPerSolarMass;
    return m * std::sqrt(m) / std::sqrt(mStarSolar);
}

// --- Stellar structure ------------------------------------------------------------------

TovStar solveTov(double K, double gamma, double rhoc) {
    // dm/dr = 4 pi r^2 eps,  dP/dr = -(eps + P)(m + 4 pi r^3 P) / (r (r - 2m)).
    auto energy = [&](double P) {
        if (P <= 0.0) return 0.0;
        const double rho0 = std::pow(P / K, 1.0 / gamma);
        return rho0 + P / (gamma - 1.0);
    };
    auto rhs = [&](double r, double m, double P, double& dm, double& dP) {
        const double eps = energy(P);
        const double p = std::max(P, 0.0);
        dm = 4.0 * kPi * r * r * eps;
        dP = -(eps + p) * (m + 4.0 * kPi * r * r * r * p) / (r * (r - 2.0 * m));
    };
    const double Pc = K * std::pow(rhoc, gamma);
    const double epsc = energy(Pc);
    double r = 1e-6, h = 1e-3;
    double m = 4.0 / 3.0 * kPi * r * r * r * epsc, P = Pc;
    for (int i = 0; i < 2000000; ++i) {
        double k1m, k1p, k2m, k2p, k3m, k3p, k4m, k4p;
        rhs(r, m, P, k1m, k1p);
        rhs(r + 0.5 * h, m + 0.5 * h * k1m, P + 0.5 * h * k1p, k2m, k2p);
        rhs(r + 0.5 * h, m + 0.5 * h * k2m, P + 0.5 * h * k2p, k3m, k3p);
        rhs(r + h, m + h * k3m, P + h * k3p, k4m, k4p);
        const double mNext = m + h / 6.0 * (k1m + 2.0 * k2m + 2.0 * k3m + k4m);
        const double pNext = P + h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
        if (pNext <= 0.0) {
            // The surface lies within this step. Near it P falls like (R - r)^2 (Gamma = 2),
            // so refine with a few shrinking steps rather than a linear guess.
            if (h < 1e-9) {
                r += h * P / (P - pNext);
                break;
            }
            h *= 0.25;
            continue;
        }
        r += h;
        m = mNext;
        P = pNext;
        // Larger steps where the profile is smooth, smaller near the surface.
        h = std::min(1e-3, std::max(1e-9, 0.02 * P / std::max(std::abs(k1p), 1e-300)));
        h = std::min(h, 1e-3);
    }
    return {m, r};
}

TovStar maxTovMass(double K, double gamma, double* centralDensity) {
    // Golden-section search for the peak of M(rho_c).
    double lo = 1e-4, hi = 2e-2;
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
    double f1 = solveTov(K, gamma, x1).mass, f2 = solveTov(K, gamma, x2).mass;
    for (int i = 0; i < 60; ++i) {
        if (f1 < f2) {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + g * (hi - lo);
            f2 = solveTov(K, gamma, x2).mass;
        } else {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - g * (hi - lo);
            f1 = solveTov(K, gamma, x1).mass;
        }
    }
    const double rho = 0.5 * (lo + hi);
    if (centralDensity) *centralDensity = rho;
    return solveTov(K, gamma, rho);
}

double neutronStarMaxMass() {
    static const double mass = maxTovMass(kNeutronStarK, kNeutronStarGamma).mass;
    return mass;
}

double neutronStarRadiusKm(double massSolar) {
    static double rhoMax = 0.0;
    static const double mMax = maxTovMass(kNeutronStarK, kNeutronStarGamma, &rhoMax).mass;
    if (massSolar > mMax) return 0.0;
    // M(rho_c) rises monotonically on the stable branch: bisect on log rho_c.
    double lo = std::log(1e-6), hi = std::log(rhoMax);
    TovStar star{};
    for (int i = 0; i < 50; ++i) {
        const double mid = 0.5 * (lo + hi);
        star = solveTov(kNeutronStarK, kNeutronStarGamma, std::exp(mid));
        (star.mass < massSolar ? lo : hi) = mid;
    }
    return star.radius * kKmPerSolarMass;
}

double whiteDwarfRadiusKm(double massSolar) {
    if (massSolar >= kChandrasekharMass) return 0.0;
    const double x = std::pow(massSolar / kChandrasekharMass, 2.0 / 3.0);
    return 0.0112 * kSunRadiusKm * std::sqrt(1.0 / x - x);
}

// --- Collapse ------------------------------------------------------------------------------

bool capSpin(Object& object) {
    const double limit = object.maxChi();
    if (object.chi <= limit) return false;
    object.chi = limit;
    return true;
}

namespace {

void becomeBlackHole(Object& o, std::string& message, const std::string& why) {
    char buf[320];
    const bool shed = o.chi >= 0.998;
    o.kind = Kind::BlackHole;
    o.name = "Black hole";
    o.radius = 0.0;
    o.pulsar = o.magnetar = false;
    o.magneticField = 0.0;
    o.surfaceTemperature = 0.0;
    if (shed) {
        std::snprintf(buf, sizeof(buf),
                      "%s It collapsed to a black hole. With J/M^2 = %.3g it first had to shed angular "
                      "momentum: a black hole cannot spin faster than a = 1, and accretion stops near 0.998.",
                      why.c_str(), o.chi);
        o.chi = 0.998;
    } else {
        std::snprintf(buf, sizeof(buf), "%s It collapsed to a black hole with the same mass and spin (a = %.3f).",
                      why.c_str(), o.chi);
    }
    message = buf;
}

}  // namespace

bool collapseIfUnstable(Object& o, std::string& message) {
    if (o.hasHorizon()) return false;
    if (o.radius < kBuchdahlRadius) {
        becomeBlackHole(o, message, "Smaller than the Buchdahl limit (R = 9/4 M), no pressure can hold a star up.");
        return true;
    }
    if (o.kind == Kind::NeutronStar && o.massSolar > neutronStarMaxMass()) {
        char why[160];
        std::snprintf(why, sizeof(why), "Above the maximum neutron-star mass (%.2f solar masses for this equation of "
                      "state), neutron degeneracy pressure fails.", neutronStarMaxMass());
        becomeBlackHole(o, message, why);
        return true;
    }
    if (o.kind == Kind::WhiteDwarf && o.massSolar > kChandrasekharMass) {
        const double rKm = neutronStarRadiusKm(o.massSolar);
        if (rKm <= 0.0) {
            becomeBlackHole(o, message, "Above the Chandrasekhar mass, and too heavy even for a neutron star.");
            return true;
        }
        // Accretion-induced collapse: electrons are squeezed into the nuclei, and what is
        // left is a neutron star. J is kept; the much smaller star spins far faster.
        o.kind = Kind::NeutronStar;
        o.name = "Neutron star";
        o.radius = rKm / (kKmPerSolarMass * o.massSolar);
        o.inertia = 0.35;
        o.surfaceTemperature = 1.0e6;
        capSpin(o);
        char buf[320];
        std::snprintf(buf, sizeof(buf),
                      "Above the Chandrasekhar mass (1.44 solar masses), electron degeneracy pressure fails: the "
                      "white dwarf collapsed into a %.1f km neutron star (accretion-induced collapse).",
                      rKm);
        message = buf;
        return true;
    }
    return false;
}

bool becomeStage(Object& o, Kind stage, double originalTemperature) {
    if (o.kind == stage || o.hasHorizon()) return false;
    o.kind = stage;
    o.pulsar = o.magnetar = false;
    switch (stage) {
        case Kind::NeutronStar:
            o.name = "Neutron star";
            o.surfaceTemperature = 1.0e6;
            o.magneticField = 1e9;
            o.inertia = 0.35;
            break;
        case Kind::WhiteDwarf:
            o.name = "White dwarf";
            o.surfaceTemperature = 20000.0;
            o.magneticField = 1e6;
            o.inertia = 0.205;
            break;
        default:
            o.name = "Star";
            o.surfaceTemperature = originalTemperature;
            o.magneticField = 0.0;
            break;
    }
    return true;
}

double displayTemperature(double kelvin) {
    return kelvin <= 30000.0 ? kelvin : 30000.0 * std::pow(kelvin / 30000.0, 0.3);
}

Kind stageForCompactness(double R, Kind original) {
    if (R <= 30.0) return Kind::NeutronStar;
    if (R <= 3.0e4) return Kind::WhiteDwarf;
    return original == Kind::NeutronStar || original == Kind::WhiteDwarf || original == Kind::BlackHole ? Kind::Star
                                                                                                       : original;
}

}  // namespace central
