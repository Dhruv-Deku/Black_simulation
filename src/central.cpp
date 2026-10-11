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
    if (km < 1e7) return printf1("%.4g km", km);
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

// --- Catalogue ------------------------------------------------------------------

const std::vector<CatalogEntry>& catalog() {
    // Spins are chi = c J / (G M^2) = I Omega / M^2 from the observed rotation: Earth
    // (1 day) ~890, Jupiter (10 h) ~860, the Sun (25 days) 0.22, a 716 Hz pulsar 0.37.
    static const std::vector<CatalogEntry> list = {
        {"Earth", "Our planet: 3e-6 solar masses, 6,371 km. A spin far above 1 (J/M^2 ~ 900) is normal for a planet; "
                  "only black holes are limited to 1.",
         makeBody(Kind::RockyPlanet, "Earth", 3.0e-6, 6371.0, 890.0, 0.33, 288.0), 32.0, false},
        {"Jupiter", "The largest planet: 0.001 solar masses, 70,000 km across, spinning every 10 hours.",
         makeBody(Kind::GasGiant, "Jupiter", 9.55e-4, 69911.0, 860.0, 0.254, 165.0), 32.0, false},
        {"Brown dwarf", "A failed star, 50 Jupiters, too light to fuse hydrogen: it glows dull red from its heat of "
                        "formation (~1,500 K).",
         makeBody(Kind::Star, "Brown dwarf", 0.05, 70000.0, 6.0, 0.205, 1500.0), 32.0, false},
        {"Red dwarf", "Proxima Centauri: 0.12 solar masses, a cool (3,000 K) red star, the most common kind.",
         makeBody(Kind::Star, "Red dwarf", 0.12, 107000.0, 0.05, 0.205, 3042.0), 32.0, false},
        {"Sun", "The Sun: 695,700 km, 5,772 K, spin 0.22. Its gravitational redshift is two parts in a million.",
         makeBody(Kind::Star, "Sun", 1.0, kSunRadiusKm, 0.22, 0.07, 5772.0), 32.0, false},
        {"Red giant", "A Sun in old age, swollen to 100 times its size and cooled to 4,000 K. Its spin is the Sun's: "
                      "J/M^2 does not change as it swells.",
         makeBody(Kind::Star, "Red giant", 1.0, 100.0 * kSunRadiusKm, 0.2, 0.1, 4000.0), 32.0, false},
        {"White dwarf", "A typical white dwarf: 0.6 solar masses in an Earth-sized ball (8,700 km), held up by "
                        "degenerate electrons. Here feeding from a disk that reaches its surface (a dwarf nova). Raise the magnetic "
                        "field to make an intermediate polar.",
         makeBody(Kind::WhiteDwarf, "White dwarf", 0.6, 8700.0, 0.004, 0.205, 15000.0, 0.0), 32.0, true},
        {"Sirius B", "The companion of Sirius: 1.02 solar masses, 5,840 km, 25,000 K. Its gravitational redshift "
                     "(80 km/s) was measured with the Hubble Space Telescope.",
         makeBody(Kind::WhiteDwarf, "Sirius B", 1.02, 5840.0, 0.002, 0.15, 25000.0), 32.0, false},
        {"Neutron star", "1.4 solar masses in 12 km (R = 5.8 M): smaller than its own ISCO. Light from the back of "
                         "the star bends round to you, and the disk stops at the ISCO, not the surface. An old star "
                         "with a weak field (1e7 G), too weak to hold the disk off.",
         makeBody(Kind::NeutronStar, "Neutron star", 1.4, 12.0, 0.1, 0.35, 1.0e6, 1e7, 30.0), 32.0, true},
        {"Pulsar", "A millisecond pulsar spinning 716 times a second (chi = 0.37), its weak field (1e8 G) lets "
                   "the disk reach close in. Beams sweep from its magnetic poles.",
         [] {
             Object o = makeBody(Kind::NeutronStar, "Pulsar", 1.4, 12.0, 0.37, 0.35, 1.0e6, 1e8, 30.0);
             o.pulsar = true;
             return o;
         }(),
         32.0, true},
        {"Magnetar", "A neutron star with a 1e15 G field, turning once every 5 s. Its crust cracks under the "
                     "magnetic stress and flares.",
         [] {
             Object o = makeBody(Kind::NeutronStar, "Magnetar", 1.4, 12.0, 1.0e-4, 0.35, 3.0e6, 1e15, 45.0);
             o.magnetar = true;
             return o;
         }(),
         32.0, false},
        {"Stellar black hole", "A 10 solar-mass black hole, like those found in X-ray binaries.",
         makeHole("Stellar black hole", 10.0, 0.7), 32.0, true},
        {"Sgr A*", "The black hole at the center of the Milky Way: 4.3 million Suns. 1 M of time is 21 s.",
         makeHole("Sgr A*", 4.3e6, 0.9), 32.0, true},
        {"M87*", "The first black hole imaged: 6.5 billion Suns. 1 M is 64 AU across; 1 M of time is 9 hours.",
         makeHole("M87*", 6.5e9, 0.9), 32.0, true},
    };
    return list;
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

double eddingtonRate(double massSolar) {
    constexpr double c = 2.99792458e10;
    return 1.26e38 * massSolar / (0.1 * c * c);
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
