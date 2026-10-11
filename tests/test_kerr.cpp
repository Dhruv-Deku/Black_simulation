// Physics tests for the Kerr geodesic code, checked against closed-form results.

#include "debris.h"
#include "gravity.h"
#include "kerr.h"
#include "tidal.h"
#include "trajectory.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

using namespace kerr;

static int g_failures = 0;
static int g_checks = 0;

static void check(bool ok, const char* what, double got, double want) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::printf("FAIL  %-60s got %.12g, want %.12g\n", what, got, want);
    }
}

static void near(double got, double want, double tol, const char* what) {
    check(std::abs(got - want) <= tol, what, got, want);
}

static void truth(bool ok, const char* what) {
    check(ok, what, ok ? 1.0 : 0.0, 1.0);
}

// Complete elliptic integral of the first kind K(m), m = k^2, via the AGM.
static double ellipticK(double m) {
    double a = 1.0, b = std::sqrt(1.0 - m);
    for (int i = 0; i < 30; ++i) {
        double an = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = an;
    }
    return kPi / (2.0 * a);
}

// --- Characteristic radii -------------------------------------------------------

static void testRadii() {
    near(horizonRadius(0.0), 2.0, 1e-14, "horizon a=0");
    near(horizonRadius(0.6), 1.8, 1e-14, "horizon a=0.6");
    near(ergosphereRadius(0.9, kPi / 2), 2.0, 1e-14, "ergosphere equator");
    near(ergosphereRadius(0.6, 0.0), horizonRadius(0.6), 1e-14, "ergosphere touches horizon at pole");

    near(iscoRadius(0.0, true), 6.0, 1e-12, "ISCO a=0");
    near(iscoRadius(0.0, false), 6.0, 1e-12, "ISCO a=0 retrograde");
    near(iscoRadius(0.5, true), 4.233003, 1e-5, "ISCO a=0.5 prograde");
    near(iscoRadius(1.0, true), 1.0, 1e-12, "ISCO a=1 prograde");
    near(iscoRadius(1.0, false), 9.0, 1e-12, "ISCO a=1 retrograde");

    near(photonOrbitRadius(0.0, true), 3.0, 1e-12, "photon orbit a=0");
    near(photonOrbitRadius(1.0, true), 1.0, 1e-12, "photon orbit a=1 prograde");
    near(photonOrbitRadius(1.0, false), 4.0, 1e-12, "photon orbit a=1 retrograde");
    near(marginallyBoundRadius(0.0, true), 4.0, 1e-12, "marginally bound a=0");
}

// --- Circular orbits ---------------------------------------------------------------

static void testCircularOrbits() {
    Constants c = circularOrbit(0.0, 6.0, true);
    near(c.E, std::sqrt(8.0 / 9.0), 1e-14, "ISCO energy a=0");
    near(c.L, std::sqrt(12.0), 1e-13, "ISCO angular momentum a=0");

    // R(r) and R'(r) vanish on a circular orbit.
    for (bool pro : {true, false}) {
        double a = 0.9, r = pro ? 5.0 : 10.0;
        Constants k = circularOrbit(a, r, pro);
        double scale = r * r * r * r;
        double h = 1e-5;
        double dR = (radialPotential(a, k, r + h) - radialPotential(a, k, r - h)) / (2 * h);
        near(radialPotential(a, k, r) / scale, 0.0, 1e-13, pro ? "R=0 prograde circular" : "R=0 retrograde circular");
        near(dR / scale, 0.0, 1e-9, pro ? "R'=0 prograde circular" : "R'=0 retrograde circular");
    }

    // Launching with the ZAMO circular speed reproduces the BPT constants.
    for (bool pro : {true, false}) {
        double a = 0.9, r = 8.0;
        double v[3] = {0.0, 0.0, circularZamoSpeed(a, r, pro)};
        Constants k;
        State s;
        launchFromZamo(a, r, kPi / 2, 0.0, v, 1.0, k, s);
        Constants want = circularOrbit(a, r, pro);
        near(k.E, want.E, 1e-12, pro ? "ZAMO launch E prograde" : "ZAMO launch E retrograde");
        near(k.L, want.L, 1e-12, pro ? "ZAMO launch L prograde" : "ZAMO launch L retrograde");
    }
}

static void testZamoRoundTrip() {
    double a = 0.8;
    double v[3] = {0.21, -0.33, 0.47};
    Constants c;
    State s;
    launchFromZamo(a, 7.5, 1.1, 0.4, v, 1.0, c, s);
    double back[3];
    zamoVelocity(a, c, s, back);
    near(back[0], v[0], 1e-12, "ZAMO round trip v_r");
    near(back[1], v[1], 1e-12, "ZAMO round trip v_theta");
    near(back[2], v[2], 1e-12, "ZAMO round trip v_phi");

    // The 4-velocity is unit timelike: g_mu_nu u^mu u^nu = -1.
    auto u = fourVelocity(a, c, s);
    double norm = dot(metricTensor(a, s.r, s.theta), u, u);
    near(norm, -1.0, 1e-12, "4-velocity normalization");
}

// --- Orbit integration -------------------------------------------------------------

static void integrate(double a, const Constants& c, State& s, double tEnd, double tol) {
    while (s.t < tEnd) {
        double h = std::min(suggestedTimeStep(a, c, s, tol), tEnd - s.t);
        stepCoordinateTime(a, c, s, h);
    }
}

static void testCircularIntegration() {
    double a = 0.7, r = 8.0;
    Constants c = circularOrbit(a, r, true);
    State s{0.0, r, kPi / 2, 0.0, 0.0, 0.0};
    double omega = circularOmega(a, r, true);
    double tEnd = 10.0 * 2.0 * kPi / omega;
    integrate(a, c, s, tEnd, 0.01);
    near(s.r, r, 1e-8, "circular orbit keeps its radius over 10 orbits");
    near(s.phi, omega * tEnd, 1e-6, "circular orbit phase matches Omega t");
    near(s.theta, kPi / 2, 1e-12, "circular orbit stays equatorial");
}

static void testConstraintConservation() {
    // A generic inclined, eccentric orbit around a rapidly spinning hole.
    double a = 0.95;
    double v[3] = {0.05, 0.25, 0.30};
    Constants c;
    State s;
    launchFromZamo(a, 12.0, 1.2, 0.0, v, 1.0, c, s);
    integrate(a, c, s, 3000.0, 0.005);
    double r4 = s.r * s.r * s.r * s.r;
    near((s.vr * s.vr - radialPotential(a, c, s.r)) / r4, 0.0, 1e-8, "radial constraint after 3000 M");
    near(s.vtheta * s.vtheta - polarPotential(a, c, s.theta), 0.0, 1e-8, "polar constraint after 3000 M");
}

static void testEquatorialOrbitSolver() {
    double a = 0.6, rp = 7.0, ra = 25.0;
    double E, L;
    truth(equatorialOrbit(a, rp, ra, true, E, L), "solver finds rp=7, ra=25 orbit");
    State s{0.0, ra, kPi / 2, 0.0, 0.0, 0.0};
    Constants c{E, L, 0.0, 1.0};
    double rMin = ra, rMax = 0.0;
    // Integrate in Mino time over a little more than one radial period.
    for (int i = 0; i < 400000 && s.t < 1500.0; ++i) {
        stepMino(a, c, s, 2e-4);
        rMin = std::min(rMin, s.r);
        rMax = std::max(rMax, s.r);
    }
    near(rMin, rp, 1e-5, "solved orbit reaches periapsis");
    near(rMax, ra, 1e-5, "solved orbit returns to apoapsis");

    double rs;
    truth(separatrixPeriapsis(0.0, 30.0, true, rs), "separatrix found for a=0");
    // Schwarzschild separatrix: p = 6 + 2e, so rp = p / (1 + e).
    // With ra fixed: rp = (6 + 2e)/(1 + e) and e = (ra - rp)/(ra + rp).
    double e = (30.0 - rs) / (30.0 + rs);
    near(2.0 * 30.0 * rs / (30.0 + rs), 6.0 + 2.0 * e, 1e-6, "separatrix p = 6 + 2e");
}

static void testPeriapsisPrecession() {
    // Schwarzschild: phi advance per radial period is
    //   4 sqrt(p / (p - 6 + 2e)) K(4e / (p - 6 + 2e)).
    double p = 20.0, e = 0.3;
    double rp = p / (1 + e), ra = p / (1 - e);
    double E, L;
    truth(equatorialOrbit(0.0, rp, ra, true, E, L), "solver finds p=20, e=0.3");
    double want = 4.0 * std::sqrt(p / (p - 6.0 + 2.0 * e)) * ellipticK(4.0 * e / (p - 6.0 + 2.0 * e));

    Constants c{E, L, 0.0, 1.0};
    State s{0.0, rp, kPi / 2, 0.0, 0.0, 0.0};
    // March in Mino time; find the next periapsis (vr crosses 0 from - to +).
    double prevVr = 0.0, prevPhi = 0.0, got = 0.0;
    bool leftPeri = false;
    for (int i = 0; i < 10000000; ++i) {
        prevVr = s.vr;
        prevPhi = s.phi;
        stepMino(0.0, c, s, 1e-4);
        if (s.vr < 0.0) leftPeri = true;
        if (leftPeri && prevVr < 0.0 && s.vr >= 0.0) {
            double f = -prevVr / (s.vr - prevVr);
            got = prevPhi + f * (s.phi - prevPhi);
            break;
        }
    }
    near(got, want, 1e-6, "exact Schwarzschild periapsis advance");
}

// --- Photons and the shadow -----------------------------------------------------------

// Sends an equatorial photon inward from far away with impact parameter b = L/E.
// Returns +1 if it escapes back to infinity, -1 if it falls into the hole.
static int photonFate(double a, double b, double* measuredB) {
    double r = 1000.0, th = kPi / 2;
    Metric m = metricAt(a, r, th);
    // Solve L/E = b for the ZAMO-frame direction.
    double vphi = b * std::sqrt(m.sigma * m.delta / m.A) / (std::sqrt(m.A / m.sigma) * (1.0 - m.omega * b));
    double v[3] = {-std::sqrt(1.0 - vphi * vphi), 0.0, vphi};
    Constants c;
    State s;
    launchFromZamo(a, r, th, 0.0, v, 0.0, c, s);
    if (measuredB) *measuredB = c.L / c.E;

    double rh = horizonRadius(a);
    for (int i = 0; i < 2000000; ++i) {
        stepMino(a, c, s, photonStepSize(a, c, s, 1e-3));
        enforceConstraints(a, c, s);  // as the GPU ray tracer does
        if (s.r < rh + 1e-3) return -1;
        if (s.r > 2000.0) return +1;
    }
    return 0;
}

static void testShadow() {
    double bc = 3.0 * std::sqrt(3.0);
    double got;
    truth(photonFate(0.0, bc * 0.999, &got) == -1, "a=0: b just inside 3 sqrt(3) is captured");
    near(got, bc * 0.999, 1e-10, "photon launch reproduces impact parameter");
    truth(photonFate(0.0, bc * 1.001, nullptr) == +1, "a=0: b just outside 3 sqrt(3) escapes");

    // Kerr a = 0.9: equatorial critical impact parameters from the photon orbits,
    // xi(r) = -(r^3 - 3r^2 + a^2 r + a^2) / (a (r - 1)).
    double a = 0.9;
    auto xi = [a](double r) { return -(r * r * r - 3 * r * r + a * a * r + a * a) / (a * (r - 1.0)); };
    double bPro = xi(photonOrbitRadius(a, true));
    double bRetro = xi(photonOrbitRadius(a, false));
    near(bPro, 2.844, 2e-3, "a=0.9 prograde critical b");
    near(bRetro, -6.834, 2e-3, "a=0.9 retrograde critical b");
    truth(photonFate(a, bPro * 0.999, nullptr) == -1, "a=0.9 prograde: inside is captured");
    truth(photonFate(a, bPro * 1.001, nullptr) == +1, "a=0.9 prograde: outside escapes");
    truth(photonFate(a, bRetro * 0.999, nullptr) == -1, "a=0.9 retrograde: inside is captured");
    truth(photonFate(a, bRetro * 1.001, nullptr) == +1, "a=0.9 retrograde: outside escapes");
}

// Traces a photon launched from (r, theta, 0) along the ZAMO-frame direction v, the
// way the shader does, and returns the direction it leaves in.
static std::array<double, 3> escapeDirection(double a, double r, double theta, const double v[3], double tol) {
    Constants c;
    State s;
    launchFromZamo(a, r, theta, 0.0, v, 0.0, c, s);
    for (int i = 0; i < 100000 && s.r < 500.0; ++i) {
        stepMino(a, c, s, photonStepSize(a, c, s, tol));
        normalizeAngles(s);
        enforceConstraints(a, c, s);
        if (s.r < horizonRadius(a) + 1e-3) return {0.0, 0.0, 0.0};
    }
    Derivatives d = derivatives(a, c, s);
    std::array<double, 3> dir = cartesianVelocity(a, s.r, s.theta, s.phi, s.vr, s.vtheta, d.dphi);
    double n = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    for (double& x : dir) x /= n;
    return dir;
}

static void testNearAxisPhotons() {
    // Rays that pass the spin axis at a distance of +-0.02 M swing through ~pi in phi
    // there. Passing on either side, they must leave in nearly the same direction
    // (frame dragging shifts both alike). A step-size rule that let a step land right
    // at the turning point once kicked one of them sideways by 15 degrees and drew a
    // bright line through the sky above the hole.
    const double a = 0.94, r = 32.0, theta = 1.48;
    for (double eps : {0.0011, 0.003}) {
        double v[2][3];
        for (int k = 0; k < 2; ++k) {
            double x = -1.0, y = -0.53, z = (k == 0 ? eps : -eps);
            double n = std::sqrt(x * x + y * y + z * z);
            v[k][0] = x / n;
            v[k][1] = y / n;
            v[k][2] = z / n;
        }
        auto plus = escapeDirection(a, r, theta, v[0], 0.04);
        auto minus = escapeDirection(a, r, theta, v[1], 0.04);
        near(plus[0], minus[0], 2e-3, "near-axis rays leave together (x)");
        near(plus[1], minus[1], 1e-2, "near-axis rays leave together (y)");
        near(plus[2], minus[2], 2e-3, "near-axis rays leave together (z)");
    }

    // Looking down from 2 degrees off the axis, the central column of pixels crosses
    // the axis almost exactly. The outgoing direction must vary smoothly with the
    // offset, from the snapped (L = 0) rays through rays that just miss the axis.
    std::array<double, 3> previous{};
    for (int k = 0; k < 6; ++k) {
        double eps = 1e-6 * std::pow(10.0, k);  // 1e-6 ... 0.1
        double v[3] = {-0.8, -0.6, eps};
        double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (double& x : v) x /= n;
        auto dir = escapeDirection(a, r, 0.035, v, 0.04);
        if (k > 0) {
            double diff = std::sqrt((dir[0] - previous[0]) * (dir[0] - previous[0]) +
                                    (dir[1] - previous[1]) * (dir[1] - previous[1]) +
                                    (dir[2] - previous[2]) * (dir[2] - previous[2]));
            // Neighbouring offsets differ by 10x, so allow a change proportional to eps.
            near(diff, 0.0, 2e-3 + 2.0 * eps, "top-down rays vary smoothly near the axis");
        }
        previous = dir;
    }
}

static void testCartesianRoundTrip() {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> uni(-20.0, 20.0);
    double worst = 0.0;
    for (int i = 0; i < 1000; ++i) {
        double x = uni(rng), y = uni(rng), z = uni(rng);
        auto [r, th, ph] = fromCartesian(0.9, x, y, z);
        auto p = toCartesian(0.9, r, th, ph);
        worst = std::max({worst, std::abs(p[0] - x), std::abs(p[1] - y), std::abs(p[2] - z)});
    }
    near(worst, 0.0, 1e-10, "Cartesian <-> spheroidal round trip");
}

// Contracts a photon (constants c, state s moving toward the observer) with an
// observer's 4-velocity: the energy that observer measures.
static double measuredEnergy(double a, const Constants& c, const State& s, const Vec4& u) {
    const Derivatives d = derivatives(a, c, s);
    const double sigma = metricAt(a, s.r, s.theta).sigma;
    const Vec4 p = {d.dt / sigma, d.dr / sigma, d.dtheta / sigma, d.dphi / sigma};
    return -dot(metricTensor(a, s.r, s.theta), p, u);
}

static void testMovingObserver() {
    const double a = 0.9, r = 7.0, th = 1.2;
    double look[3] = {-0.6, 0.48, 0.64};  // unit

    // At rest relative to the ZAMO there is no aberration.
    double rest[3] = {0.0, 0.0, 0.0};
    Constants c0;
    State s0;
    double g0 = observedPhoton(a, r, th, 0.0, rest, look, c0, s0);

    // Moving: the measured energy equals what the boost formula predicts, and the
    // photon is still null (vr^2 = R, vtheta^2 = Theta).
    double vObs[3] = {-0.3, 0.2, 0.5};
    Constants c;
    State s;
    double g = observedPhoton(a, r, th, 0.0, vObs, look, c, s);
    Constants oc;
    State os;
    launchFromZamo(a, r, th, 0.0, vObs, 1.0, oc, os);
    near(measuredEnergy(a, c, s, fourVelocity(a, oc, os)), g, 1e-10, "moving observer measures the predicted energy");
    near(s.vr * s.vr - radialPotential(a, c, s.r), 0.0, 1e-9, "aberrated photon is null (radial)");
    near(s.vtheta * s.vtheta - polarPotential(a, c, s.theta), 0.0, 1e-9, "aberrated photon is null (polar)");

    // Flying straight at what you look at blueshifts it by sqrt((1 + v) / (1 - v)).
    double ahead[3] = {-1.0, 0.0, 0.0};
    double forward[3] = {-0.6, 0.0, 0.0};
    Constants c1, c2;
    State s1, s2;
    double still = observedPhoton(a, r, th, 0.0, rest, ahead, c1, s1);
    double moving = observedPhoton(a, r, th, 0.0, forward, ahead, c2, s2);
    near(moving / still, 2.0, 1e-12, "head-on Doppler factor sqrt(1.6/0.4) = 2");
    near(g0, observedPhoton(a, r, th, 0.0, rest, look, c1, s1), 0.0, "observer at rest is deterministic");
}

static void testDebris() {
    std::mt19937 rng(3);
    // A cloud at rest (a = 0): each piece has E = sqrt(1 - 2/r_i), so a ball of
    // radius R at distance r spans about 2R / (r^2 sqrt(1 - 2/r)) in energy. This
    // spread is what pulls a disrupted star into a stream.
    double still[3] = {0.0, 0.0, 0.0};
    auto parts = debris::cloud(0.0, 20.0, kPi / 2, 0.0, still, 1.0, false, 4000, 10000.0f, rng);
    truth(parts.size() == 4000, "cloud has the requested particle count");
    double worst = 0.0, eMin = 2.0, eMax = 0.0;
    for (const ParticleInit& p : parts) {
        worst = std::max(worst, std::abs(p.E - std::sqrt(1.0 - 2.0 / p.r)));
        eMin = std::min(eMin, static_cast<double>(p.E));
        eMax = std::max(eMax, static_cast<double>(p.E));
    }
    near(worst, 0.0, 2e-6, "static cloud: E = sqrt(1 - 2/r) for every piece");
    double predicted = 2.0 / (400.0 * std::sqrt(1.0 - 2.0 / 20.0));
    near((eMax - eMin) / predicted, 1.0, 0.05, "tidal energy spread 2R / (r^2 sqrt(1 - 2/r))");

    // A tiny cloud moving like an orbiting body keeps that body's orbit.
    double v[3] = {0.1, -0.05, 0.3};
    auto tiny = debris::cloud(0.7, 15.0, 1.3, 0.5, v, 1e-5, false, 50, 5000.0f, rng);
    Constants c;
    State s;
    launchFromZamo(0.7, 15.0, 1.3, 0.5, v, 1.0, c, s);
    double dE = 0.0, dL = 0.0;
    for (const ParticleInit& p : tiny) {
        dE = std::max(dE, std::abs(p.E - c.E));
        dL = std::max(dL, std::abs(p.L - c.L));
    }
    near(dE, 0.0, 1e-5, "tiny moving cloud: same energy as its center");
    near(dL, 0.0, 1e-4, "tiny moving cloud: same angular momentum as its center");

    // A ring tilted by 30 degrees around a non-spinning hole: each orbit's
    // inclination, atan(sqrt(Q) / L), is 30 degrees.
    auto ring = debris::ring(0.0, 10.0, 12.0, 30.0, 2000, 4000.0f, rng);
    double meanTilt = 0.0;
    for (const ParticleInit& p : ring) meanTilt += std::atan2(std::sqrt(std::max(p.Q, 0.0f)), p.L);
    meanTilt /= static_cast<double>(ring.size());
    near(meanTilt * 180.0 / kPi, 30.0, 0.5, "ring orbits are tilted by 30 degrees");
}

static void testTrajectoryClassification() {
    using trajectory::Fate;
    // Schwarzschild effective potential: V = sqrt((1 - 2/r)(1 + L^2 / r^2)).
    double worst = 0.0;
    for (double r = 2.5; r < 50.0; r *= 1.3) {
        worst = std::max(worst, std::abs(trajectory::effectivePotential(0.0, 4.0, 0.0, r) -
                                         std::sqrt((1.0 - 2.0 / r) * (1.0 + 16.0 / (r * r)))));
    }
    near(worst, 0.0, 1e-12, "effective potential matches Schwarzschild");

    // Each case: classification, and the numerical prediction agrees with it.
    struct Case {
        const char* name;
        double a, r;
        double v[3];
        Fate want;
    };
    auto launch = [](double a, double r, const double v[3], Constants& c, State& s) {
        launchFromZamo(a, r, kPi / 2, 0.0, v, 1.0, c, s);
    };
    const Case cases[] = {
        {"plunge (slow orbit)", 0.9, 14.0, {0.0, 0.0, 0.45 * circularZamoSpeed(0.9, 14.0, true)}, Fate::Plunges},
        {"escape (hyperbolic flyby)", 0.5, 40.0, {-0.27, 0.0, 0.13}, Fate::Escapes},
        {"thrown outward, falls back", 0.0, 10.0, {0.3, 0.0, 0.0}, Fate::Plunges},
        {"thrown outward, escapes", 0.0, 10.0, {0.6, 0.0, 0.0}, Fate::Escapes},
        {"tilted bound orbit", 0.9, 12.0, {0.05, 0.2, 0.25}, Fate::Bound},
    };
    for (const Case& k : cases) {
        Constants c;
        State s;
        launch(k.a, k.r, k.v, c, s);
        auto an = trajectory::classify(k.a, c, s);
        check(an.fate == k.want, k.name, static_cast<double>(an.fate), static_cast<double>(k.want));
        auto path = trajectory::predict(k.a, c, s, 20000.0, 1.0);
        bool consistent = k.want == Fate::Bound ? path.fate == Fate::Bound : path.fate == k.want;
        check(consistent, "prediction agrees with classification", static_cast<double>(path.fate),
              static_cast<double>(k.want));
    }

    // Bound orbit with known turning points.
    double E, L;
    equatorialOrbit(0.6, 7.0, 25.0, true, E, L);
    Constants c{E, L, 0.0, 1.0};
    State s{0.0, 25.0, kPi / 2, 0.0, 0.0, 0.0};
    auto an = trajectory::classify(0.6, c, s);
    truth(an.fate == Fate::Bound, "rp=7, ra=25 orbit is bound");
    near(an.rInner, 7.0, 1e-6, "inner turning point found");
    near(an.rOuter, 25.0, 1e-5, "outer turning point found (it starts there)");

    // Circular orbits: stable outside the ISCO, unstable inside it.
    Constants stable = circularOrbit(0.0, 10.0, true), unstable = circularOrbit(0.0, 4.5, true);
    State s10{0.0, 10.0, kPi / 2, 0.0, 0.0, 0.0}, s45{0.0, 4.5, kPi / 2, 0.0, 0.0, 0.0};
    auto a10 = trajectory::classify(0.0, stable, s10);
    auto a45 = trajectory::classify(0.0, unstable, s45);
    truth(a10.fate == Fate::Circular && a10.stableCircular, "r = 10 circular orbit is stable");
    truth(a45.fate == Fate::Circular && !a45.stableCircular, "r = 4.5 circular orbit is unstable");

    // Measured periapsis precession matches the exact Schwarzschild value.
    double p = 20.0, e = 0.3;
    equatorialOrbit(0.0, p / (1 + e), p / (1 - e), true, E, L);
    Constants cp{E, L, 0.0, 1.0};
    State sp{0.0, p / (1 + e), kPi / 2, 0.0, 0.0, 0.0};
    auto path = trajectory::predict(0.0, cp, sp, 5000.0, 0.5);
    double exact = (4.0 * std::sqrt(p / (p - 6.0 + 2.0 * e)) * ellipticK(4.0 * e / (p - 6.0 + 2.0 * e)) - 2.0 * kPi) *
                   180.0 / kPi;
    near(path.precessionDeg, exact, 0.5, "predicted precession per orbit (degrees)");
}

static Frame zamoAt(double a, double r, double theta) {
    return observerFrame(a, r, theta, zamoFourVelocity(a, r, theta));
}

static void testRayTrace() {
    trajectory::RayScene scene;
    scene.disk = false;
    double atHole[3] = {-1.0, 0.0, 0.0}, sideways[3] = {0.0, 0.0, 1.0};
    auto hole = trajectory::traceRay(0.9, 30.0, kPi / 2, 0.0, zamoAt(0.9, 30.0, kPi / 2), atHole, scene);
    auto away = trajectory::traceRay(0.9, 30.0, kPi / 2, 0.0, zamoAt(0.9, 30.0, kPi / 2), sideways, scene);
    truth(hole.end == trajectory::RayEnd::Horizon, "looking at the hole: the ray ends in the horizon");
    truth(away.end == trajectory::RayEnd::Sky, "looking sideways: the ray reaches the sky");

    // Light-travel time along an outgoing radial ray (a = 0): in Kerr-Schild time
    // dt/dr = (r + 2) / (r - 2), so the delay is [r + 4 ln(r - 2)] between the two
    // ends. It diverges at the horizon, which is why distant observers see infalling
    // things freeze there. A sphere on the line of sight at r ~ 10 gives a clean end.
    trajectory::RayScene withBody = scene;
    withBody.spheres.push_back({{10.0, 0.0, 0.0}, 0.5});
    auto radial = trajectory::traceRay(0.0, 30.0, kPi / 2, 0.0, zamoAt(0.0, 30.0, kPi / 2), atHole, withBody);
    auto tortoise = [](double x) { return x + 4.0 * std::log(x - 2.0); };
    truth(radial.end == trajectory::RayEnd::Body, "ray along the line of sight hits the sphere");
    near(radial.rEnd, 10.5, 1e-2, "it hits the sphere's near side");
    near(radial.delay, tortoise(30.0) - tortoise(radial.rEnd), 1e-3, "light-travel time along a radial ray");

    // Just outside the critical impact parameter the light circles the hole first.
    const double r = 1000.0, bc = 3.0 * std::sqrt(3.0);
    auto lookFor = [&](double b, double out[3]) {
        double sa = b * std::sqrt(1.0 - 2.0 / r) / r;
        out[0] = -std::sqrt(1.0 - sa * sa);
        out[1] = 0.0;
        out[2] = sa;
    };
    double grazing[3], wide[3];
    lookFor(bc * (1.0 + 1e-6), grazing);
    lookFor(2.0 * bc, wide);
    auto g1 = trajectory::traceRay(0.0, r, kPi / 2, 0.0, zamoAt(0.0, r, kPi / 2), grazing, scene);
    auto g2 = trajectory::traceRay(0.0, r, kPi / 2, 0.0, zamoAt(0.0, r, kPi / 2), wide, scene);
    truth(g1.end == trajectory::RayEnd::Sky && g1.turns > 1.5, "grazing ray circles the hole before escaping");
    near(g1.rMin, 3.0, 0.01, "grazing ray skims the photon sphere at r = 3");
    truth(g2.turns < 0.7, "wide ray barely bends");

    // With a disk: looking just below the hole from above the disk lands on it,
    // and the measured frequency shift is sensible.
    scene.disk = true;
    scene.diskInner = iscoRadius(0.9, true);
    scene.diskOuter = 20.0;
    double down[3] = {-0.95, 0.3122, 0.0};
    auto disk = trajectory::traceRay(0.9, 30.0, 1.35, 0.0, zamoAt(0.9, 30.0, 1.35), down, scene);
    truth(disk.end == trajectory::RayEnd::Disk, "looking down at the disk: the ray ends on the disk");
    truth(disk.g > 0.2 && disk.g < 3.0, "disk frequency shift is plausible");
}


// --- Kerr-Schild coordinates and observers -------------------------------------------

static void invert4(const double m[4][4], double out[4][4]) {
    double a[4][8];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            a[i][j] = m[i][j];
            a[i][j + 4] = i == j ? 1.0 : 0.0;
        }
    }
    for (int c = 0; c < 4; ++c) {
        int p = c;
        for (int i = c + 1; i < 4; ++i) {
            if (std::abs(a[i][c]) > std::abs(a[p][c])) p = i;
        }
        for (int j = 0; j < 8; ++j) std::swap(a[c][j], a[p][j]);
        const double d = a[c][c];
        for (int j = 0; j < 8; ++j) a[c][j] /= d;
        for (int i = 0; i < 4; ++i) {
            if (i == c) continue;
            const double f = a[i][c];
            for (int j = 0; j < 8; ++j) a[i][j] -= f * a[c][j];
        }
    }
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) out[i][j] = a[i][j + 4];
    }
}

static void testKerrSchildMetric() {
    // The inverse metric has the known closed form, also inside the horizon.
    double worst = 0.0;
    for (double r : {12.0, 3.0, 1.5, 1.2, 0.8}) {
        const double a = 0.9, th = 1.1;
        const MetricTensor g = metricTensor(a, r, th);
        double inv[4][4];
        invert4(g.g, inv);
        const double sigma = r * r + a * a * std::cos(th) * std::cos(th), delta = r * r - 2 * r + a * a;
        const double s2 = std::sin(th) * std::sin(th);
        const double want[4][4] = {{-(1 + 2 * r / sigma), 2 * r / sigma, 0, 0},
                                   {2 * r / sigma, delta / sigma, 0, a / sigma},
                                   {0, 0, 1 / sigma, 0},
                                   {0, a / sigma, 0, 1 / (sigma * s2)}};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) worst = std::max(worst, std::abs(inv[i][j] - want[i][j]));
        }
    }
    near(worst, 0.0, 1e-12, "Kerr-Schild inverse metric");

    // In the Cartesian chart the spatial metric is flat plus one null term:
    // gamma_ij = delta_ij + f l_i l_j, f = 2 r^3 / (r^4 + a^2 z^2).
    worst = 0.0;
    for (double r : {8.0, 2.5, 1.2}) {
        const double a = 0.8, th = 0.9, ph = 0.7, h = 1e-6;
        // Jacobian d(x, y, z) / d(r, theta, phi) by central differences.
        double J[3][3];
        const double base[3] = {r, th, ph};
        for (int k = 0; k < 3; ++k) {
            double up[3] = {base[0], base[1], base[2]}, dn[3] = {base[0], base[1], base[2]};
            up[k] += h;
            dn[k] -= h;
            auto pu = toCartesian(a, up[0], up[1], up[2]);
            auto pd = toCartesian(a, dn[0], dn[1], dn[2]);
            for (int i = 0; i < 3; ++i) J[i][k] = (pu[i] - pd[i]) / (2 * h);
        }
        // gamma in (r, theta, phi) from the metric, mapped to Cartesian: J^-T gamma J^-1.
        const MetricTensor g = metricTensor(a, r, th);
        double m4[4][4] = {}, inv4[4][4];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) m4[i][j] = J[i][j];
        }
        m4[3][3] = 1.0;
        invert4(m4, inv4);
        auto p = toCartesian(a, r, th, ph);
        const double f = 2 * r * r * r / (r * r * r * r + a * a * p[2] * p[2]);
        const double l[3] = {(r * p[0] + a * p[1]) / (r * r + a * a), (r * p[1] - a * p[0]) / (r * r + a * a),
                             p[2] / r};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double sum = 0.0;
                for (int k = 0; k < 3; ++k) {
                    for (int m = 0; m < 3; ++m) sum += inv4[k][i] * g.g[k + 1][m + 1] * inv4[m][j];
                }
                worst = std::max(worst, std::abs(sum - ((i == j ? 1.0 : 0.0) + f * l[i] * l[j])));
            }
        }
    }
    near(worst, 0.0, 1e-6, "Cartesian chart: spatial metric is flat + f l l");

    // Outside the horizon the rates differ from Boyer-Lindquist ones by exactly the
    // coordinate shifts: dphi = dphi_BL + a vr / Delta, dt = dt_BL + (...) vr.
    const double a = 0.7;
    double v[3] = {-0.2, 0.1, 0.3};
    Constants c;
    State s;
    launchFromZamo(a, 6.0, 1.2, 0.0, v, 1.0, c, s);
    const Derivatives d = derivatives(a, c, s);
    const double delta = s.r * s.r - 2 * s.r + a * a, rra = s.r * s.r + a * a;
    const double P = c.E * rra - a * c.L, sn2 = std::sin(s.theta) * std::sin(s.theta);
    const double dphiBL = a * P / delta - a * c.E + c.L / sn2;
    const double dtBL = rra * P / delta - a * (a * c.E * sn2 - c.L);
    near(d.dphi - dphiBL, a * s.vr / delta, 1e-10, "dphi = dphi_BL + a vr / Delta");
    near(d.dt - dtBL, (rra / delta - 1.0) * s.vr, 1e-10, "dt = dt_BL + ((r^2 + a^2) / Delta - 1) vr");
    // And the shifts are the integrals of those factors.
    const double h = 1e-5;
    near((azimuthShift(a, 6.0 + h) - azimuthShift(a, 6.0 - h)) / (2 * h), a / delta, 1e-8,
         "azimuth shift derivative");
    near((timeShift(a, 6.0 + h) - timeShift(a, 6.0 - h)) / (2 * h), rra / delta - 1.0, 1e-8,
         "time shift derivative");
}

static void testFallThroughHorizon() {
    // Radial free fall from rest at r0 = 10 (a = 0): proper time follows the cycloid
    // r = (r0 / 2)(1 + cos eta), tau = sqrt(r0^3 / 8)(eta + sin eta), straight through
    // the horizon at r = 2 and on toward the singularity.
    const double r0 = 10.0;
    double still[3] = {0.0, 0.0, 0.0};
    Constants c;
    State s;
    launchFromZamo(0.0, r0, kPi / 2, 0.0, still, 1.0, c, s);
    auto cycloidTau = [&](double r) {
        const double eta = std::acos(2.0 * r / r0 - 1.0);
        return std::sqrt(r0 * r0 * r0 / 8.0) * (eta + std::sin(eta));
    };
    double tau = 0.0, tAtHorizon = -1.0, tauAtHorizon = 0.0, worst = 0.0;
    for (int i = 0; i < 4000000 && s.r > 0.5; ++i) {
        // Proper time: d tau = Sigma d(Mino) = r^2 d(Mino) on the equator (trapezoid).
        const State before = s;
        const double h = 2e-4;
        stepMino(0.0, c, s, h);
        tau += 0.5 * (before.r * before.r + s.r * s.r) * h;
        if (before.r >= 2.0 && s.r < 2.0) {
            const double f = (before.r - 2.0) / (before.r - s.r);
            tAtHorizon = before.t + f * (s.t - before.t);
            tauAtHorizon = tau - (1.0 - f) * 0.5 * (before.r * before.r + s.r * s.r) * h;
        }
        if (i % 1000 == 0 && s.r < r0 * 0.999) worst = std::max(worst, std::abs(tau - cycloidTau(s.r)));
    }
    near(worst, 0.0, 1e-4, "free fall: proper time matches the cycloid, through the horizon");
    truth(std::isfinite(tAtHorizon) && tAtHorizon > 0.0 && tAtHorizon < 100.0, "free fall: horizon crossed at finite t");
    near(tauAtHorizon, cycloidTau(2.0), 1e-3, "free fall: proper time at the horizon");

    // A spinning hole: coordinate-time steps carry an infalling body through r+.
    const double a = 0.9;
    double inward[3] = {-0.3, 0.05, 0.1};
    launchFromZamo(a, 6.0, 1.3, 0.0, inward, 1.0, c, s);
    bool crossed = false;
    for (int i = 0; i < 100000 && s.r > interiorLimit(a); ++i) {
        stepCoordinateTime(a, c, s, suggestedTimeStep(a, c, s, 0.005));
        normalizeAngles(s);
        if (s.r < horizonRadius(a)) crossed = true;
        if (!std::isfinite(s.r) || !std::isfinite(s.t) || !std::isfinite(s.phi)) break;
    }
    truth(crossed && std::isfinite(s.t) && std::isfinite(s.phi), "a=0.9: body crosses the horizon, t and phi finite");
    near(s.vr * s.vr - radialPotential(a, c, s.r), 0.0, 1e-7, "a=0.9: radial constraint holds inside the horizon");
}

static void testObserverFrames() {
    const double a = 0.9;
    auto orthonormality = [&](double r, double th, const Frame& f) {
        const MetricTensor g = metricTensor(a, r, th);
        double worst = 0.0;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                const double want = i == j ? (i == 0 ? -1.0 : 1.0) : 0.0;
                worst = std::max(worst, std::abs(dot(g, f.e[i], f.e[j]) - want));
            }
        }
        return worst;
    };
    near(orthonormality(6.0, 1.1, zamoAt(a, 6.0, 1.1)), 0.0, 1e-12, "ZAMO frame is orthonormal");

    // A moving observer, and one falling inside the horizon.
    double v[3] = {-0.4, 0.2, 0.35};
    Constants c;
    State s;
    launchFromZamo(a, 5.0, 1.0, 0.3, v, 1.0, c, s);
    near(orthonormality(s.r, s.theta, observerFrame(a, s.r, s.theta, fourVelocity(a, c, s))), 0.0, 1e-12,
         "moving observer's frame is orthonormal");
    for (int i = 0; i < 100000 && s.r > 1.2; ++i) {
        stepCoordinateTime(a, c, s, suggestedTimeStep(a, c, s, 0.005));
        normalizeAngles(s);
    }
    truth(s.r < horizonRadius(a), "observer is inside the horizon");
    const Frame inside = observerFrame(a, s.r, s.theta, fourVelocity(a, c, s));
    near(orthonormality(s.r, s.theta, inside), 0.0, 1e-10, "infalling observer's frame is orthonormal inside");

    // Photons built from a frame are null; from the ZAMO frame they match launchFromZamo.
    double look[3] = {-0.48, 0.6, 0.64};
    Constants pc;
    State ps;
    photonFromFrame(a, s.r, s.theta, s.phi, inside, look, pc, ps);
    near(ps.vr * ps.vr - radialPotential(a, pc, ps.r), 0.0, 1e-9, "photon seen inside the horizon is null (radial)");
    near(ps.vtheta * ps.vtheta - polarPotential(a, pc, ps.theta), 0.0, 1e-9,
         "photon seen inside the horizon is null (polar)");

    photonFromFrame(a, 7.0, 1.2, 0.0, zamoAt(a, 7.0, 1.2), look, pc, ps);
    const double dir[3] = {-look[0], -look[1], -look[2]};
    Constants zc;
    State zs;
    launchFromZamo(a, 7.0, 1.2, 0.0, dir, 0.0, zc, zs);
    near(pc.L / pc.E, zc.L, 1e-12, "ZAMO-frame photon: same L/E as launchFromZamo");
    near(pc.Q / (pc.E * pc.E), zc.Q, 1e-12, "ZAMO-frame photon: same Q/E^2 as launchFromZamo");
    near(ps.vr / pc.E, zs.vr, 1e-12, "ZAMO-frame photon: same radial velocity");

    // The ray tracer's frequency shift (trace.frag, frequencyRatio) matches -p.u from the
    // metric, outside and inside the horizon. It works from the backward-traced
    // velocities and regularizes the horizon's 0/0.
    {
        auto shaderNu = [&](const Constants& pc, const State& ps, const Vec4& u) {
            const double r = ps.r, delta = r * r - 2 * r + a * a, P = pc.E * (r * r + a * a) - a * pc.L;
            const double x = pc.L - a * pc.E, K = pc.Q + x * x;
            const double vrBack = -ps.vr, vthBack = -ps.vtheta;
            const double F = inflowFactor(P, -vrBack, delta, K);
            return pc.E * u[0] - (F - pc.E) * u[1] + vthBack * u[2] - pc.L * u[3];
        };
        double worst = 0.0;
        for (int k = 0; k < 2; ++k) {
            State bs = k == 0 ? State{0.0, 6.0, 1.2, 0.4, 0.0, 0.0} : s;  // outside; then the observer inside
            Constants bc = c;
            if (k == 0) {
                double bv[3] = {0.3, -0.2, 0.4};
                launchFromZamo(a, bs.r, bs.theta, bs.phi, bv, 1.0, bc, bs);
            }
            const Vec4 u = fourVelocity(a, bc, bs);
            const Frame f = observerFrame(a, bs.r, bs.theta, u);
            // A photon some other observer sees at the same place, so that p.u is not trivially 1.
            double vOther[3] = {k == 0 ? -0.5 : 0.0, 0.1, 0.2};
            Constants oc;
            State os;
            Frame other = f;
            if (k == 0) {
                launchFromZamo(a, bs.r, bs.theta, bs.phi, vOther, 1.0, oc, os);
                other = observerFrame(a, bs.r, bs.theta, fourVelocity(a, oc, os));
            }
            for (double lx : {-0.6, 0.0, 0.7}) {
                double lk[3] = {lx, 0.5, std::sqrt(std::max(0.0, 1.0 - lx * lx - 0.25))};
                Constants qc;
                State qs;
                photonFromFrame(a, bs.r, bs.theta, bs.phi, other, lk, qc, qs);
                worst = std::max(worst, std::abs(shaderNu(qc, qs, u) - measuredEnergy(a, qc, qs, u)));
            }
        }
        near(worst, 0.0, 1e-9, "ray tracer's frequency shift matches -p.u (also inside the horizon)");
    }

    // Looking up from inside the horizon you still see the sky; looking down, black.
    trajectory::RayScene scene;
    scene.disk = false;
    double up[3] = {1.0, 0.0, 0.0}, down[3] = {-1.0, 0.0, 0.0};
    auto skyward = trajectory::traceRay(a, s.r, s.theta, s.phi, inside, up, scene);
    auto downward = trajectory::traceRay(a, s.r, s.theta, s.phi, inside, down, scene);
    truth(skyward.end == trajectory::RayEnd::Sky, "inside the horizon: looking up, the sky");
    truth(downward.end == trajectory::RayEnd::Horizon, "inside the horizon: looking down, darkness");
}

// --- Tidal deformation (affine star model) ---------------------------------------------

static void testTidalModel() {
    // Index symbols: 2/3 each for a sphere; they always sum to 2; and the closed form
    // for a prolate spheroid, A1 = (1 - e^2)/e^3 (ln((1 + e)/(1 - e)) - 2e).
    auto sphere = tidal::indexSymbols({1.0, 1.0, 1.0});
    near(sphere[0], 2.0 / 3.0, 1e-12, "index symbols of a sphere");
    auto any = tidal::indexSymbols({0.6, 1.3, 2.4});
    near(any[0] + any[1] + any[2], 2.0, 1e-10, "index symbols sum to 2");
    const double e = 0.8;
    auto prolate = tidal::indexSymbols({1.0, std::sqrt(1.0 - e * e), std::sqrt(1.0 - e * e)});
    near(prolate[0], (1.0 - e * e) / (e * e * e) * (std::log((1.0 + e) / (1.0 - e)) - 2.0 * e), 1e-10,
         "index symbol of a prolate spheroid");

    // Far from any tide the star stays round; nudged, it oscillates without growing.
    const tidal::Mat3 none{};
    tidal::Star star;
    tidal::advance(star, none, 1.0, 200.0);
    near(tidal::stretch(star.F), 1.0, 1e-9, "an isolated star stays in equilibrium");
    star.F[0][0] = 1.05;
    double most = 0.0;
    for (int i = 0; i < 200; ++i) {
        tidal::advance(star, none, 1.0, 0.5);
        most = std::max(most, tidal::stretch(star.F));
    }
    truth(most > 1.04 && most < 1.08, "a nudged star oscillates without growing");

    // The tidal tensor is traceless and stretches along the line to the hole.
    double v[3] = {0.0, 0.0, circularZamoSpeed(0.0, 20.0, true)};
    Constants c;
    State s;
    launchFromZamo(0.0, 20.0, kPi / 2, 0.0, v, 1.0, c, s);
    const tidal::Mat3 C = tidal::tidalTensor(0.0, c, s);
    near(C[0][0] + C[1][1] + C[2][2], 0.0, 1e-15, "tidal tensor is traceless");
    // Circular orbit at r = 20: radial eigenvalue (2r - 3) / (r^3 (r - 3)).
    near(C[0][0] * 8000.0, (2.0 * 20.0 - 3.0) / (20.0 - 3.0), 1e-9, "relativistic radial tide on a circular orbit");

    // A star on the tidal-disruption preset's orbit (beta = 2.3) is torn apart; the
    // grazing preset's (beta = 0.86) survives, squeezed and wobbling.
    auto fly = [](double a, double rp, double ra, double rt, double& peak, double& after) {
        double E, L;
        equatorialOrbit(a, rp, ra, true, E, L);
        Constants k{E, L, 0.0, 1.0};
        State st{0.0, ra, kPi / 2, 0.0, 0.0, 0.0};
        tidal::Star body;
        peak = 1.0;
        bool pastPeriapsis = false;
        double prevVr = 0.0;
        for (int i = 0; i < 200000; ++i) {
            const double rate = metricAt(a, st.r, st.theta).sigma / derivatives(a, k, st).dt;
            const double h = suggestedTimeStep(a, k, st, 0.01);
            stepCoordinateTime(a, k, st, h);
            tidal::advance(body, tidal::tidalTensor(a, k, st), tidal::frequencySquared(rt), h * rate);
            peak = std::max(peak, tidal::stretch(body.F));
            if (prevVr < 0.0 && st.vr >= 0.0) pastPeriapsis = true;
            prevVr = st.vr;
            if ((pastPeriapsis && st.r > 0.9 * ra) || peak > 4.0) break;
        }
        after = tidal::stretch(body.F);
    };
    double peak, after;
    fly(0.94, 7.0, 40.0, 16.0, peak, after);
    truth(peak > 2.5, "a star plunging deep inside its tidal radius is disrupted");
    fly(0.94, 14.0, 45.0, 12.0, peak, after);
    truth(peak > 1.3 && peak < 2.0 && after < 2.0, "a star grazing its tidal radius survives, deformed");
}

// --- Gravity between bodies ----------------------------------------------------------

static void testMutualGravity() {
    // Two stars of 0.01 M each, 2 M apart, far from the hole on a mutual circular
    // orbit: after one Kepler period, 2 pi sqrt(d^3 / m), they are back where they started.
    const double a = 0.5, r = 1000.0, d = 2.0, m = 0.01;
    const double w = std::sqrt(2.0 * m / d);
    std::vector<gravity::Member> pair(2);
    for (int k = 0; k < 2; ++k) {
        const double side = k == 0 ? -1.0 : 1.0;
        double v[3] = {0.0, 0.0, -side * 0.5 * w};
        launchFromZamo(a, r + side * 0.5 * d, kPi / 2, 0.0, v, 1.0, pair[k].c, pair[k].s);
        pair[k].mass = m;
        pair[k].softening = 0.01;
    }
    const auto start0 = toCartesian(a, pair[0].s.r, pair[0].s.theta, pair[0].s.phi);
    const auto start1 = toCartesian(a, pair[1].s.r, pair[1].s.theta, pair[1].s.phi);
    const double period = 2.0 * kPi * std::sqrt(d * d * d / (2.0 * m));
    auto paths = gravity::predict(a, pair, period, 1e-3, 1000000);
    // The separation stays d all the way round, and after one period the pair is back
    // where it started.
    double worstSeparation = 0.0;
    for (size_t q = 0; q < paths[0].points.size(); ++q) {
        double d2 = 0.0;
        for (int k = 0; k < 3; ++k) d2 += (paths[1].points[q][k] - paths[0].points[q][k]) * (paths[1].points[q][k] - paths[0].points[q][k]);
        worstSeparation = std::max(worstSeparation, std::abs(std::sqrt(d2) - d));
    }
    near(worstSeparation / d, 0.0, 0.01, "binary keeps its separation (circular mutual orbit)");
    const auto& end0 = paths[0].points.back();
    const auto& end1 = paths[1].points.back();
    const double angle = std::atan2(end1[1] - end0[1], end1[0] - end0[0]) - std::atan2(start1[1] - start0[1], start1[0] - start0[0]);
    // To about 1%: our circular setup is Newtonian, the orbit relativistic.
    near(angle / (2.0 * kPi), 0.0, 0.012, "binary returns after one Kepler period");
    const auto& mid0 = paths[0].points[paths[0].points.size() / 2];
    const auto& mid1 = paths[1].points[paths[1].points.size() / 2];
    double half = 0.0;
    for (int k = 0; k < 3; ++k) half += (mid1[k] - mid0[k]) * (start1[k] - start0[k]);
    truth(half < -0.8 * d * d, "binary has swapped sides halfway round");

    // A moon at a third of its planet's Hill radius, the pair orbiting the hole at
    // r = 22 M at 0.22c, stays bound. (It spun out within a few orbits while the pull
    // ignored the relativity of simultaneity.)
    {
        const double spin = 0.94, rOrbit = 22.0, mp = 0.01, dm = 1.1;
        std::vector<gravity::Member> sys(2);
        double vp[3] = {0.0, 0.0, circularZamoSpeed(spin, rOrbit, true)};
        launchFromZamo(spin, rOrbit, kPi / 2, 0.0, vp, 1.0, sys[0].c, sys[0].s);
        sys[0].mass = mp;
        sys[0].softening = 0.06;
        double u[3] = {0.0, 0.0, std::sqrt(mp / dm)}, vm[3];
        gravity::addVelocities(vp, u, vm);
        launchFromZamo(spin, rOrbit + dm, kPi / 2, 0.0, vm, 1.0, sys[1].c, sys[1].s);
        sys[1].mass = 1e-6;
        sys[1].softening = 0.02;
        auto orbit = gravity::predict(spin, sys, 1500.0, 1e-3, 1000000);
        double widest = 0.0;
        const size_t n = std::min(orbit[0].points.size(), orbit[1].points.size());
        for (size_t q = 0; q < n; ++q) {
            double d2 = 0.0;
            for (int k = 0; k < 3; ++k) d2 += (orbit[1].points[q][k] - orbit[0].points[q][k]) * (orbit[1].points[q][k] - orbit[0].points[q][k]);
            widest = std::max(widest, std::sqrt(d2));
        }
        truth(widest < 1.3 * dm, "a moon of a planet moving at 0.22c stays bound for 1500 M");
    }
}


int main() {
    testTidalModel();
    testMutualGravity();
    testKerrSchildMetric();
    testFallThroughHorizon();
    testObserverFrames();
    testTrajectoryClassification();
    testRayTrace();
    testCartesianRoundTrip();
    testMovingObserver();
    testDebris();
    testNearAxisPhotons();
    testRadii();
    testCircularOrbits();
    testZamoRoundTrip();
    testCircularIntegration();
    testConstraintConservation();
    testEquatorialOrbitSolver();
    testPeriapsisPrecession();
    testShadow();
    std::printf("%d/%d checks passed\n", g_checks - g_failures, g_checks);
    return g_failures == 0 ? 0 : 1;
}
