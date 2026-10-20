// Physics tests for the Kerr geodesic code, checked against closed-form results.

#include "bodies.h"
#include "central.h"
#include "debris.h"
#include "eht.h"
#include "gravity.h"
#include "json.h"
#include "kerr.h"
#include "keyframes.h"
#include "merger.h"
#include "penrose.h"
#include "scene_file.h"
#include "surface.h"
#include "tde.h"
#include "tidal.h"
#include "timeline.h"
#include "trajectory.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <thread>

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


// --- Central objects: neutron stars, white dwarfs, stars, planets ---------------------

static void testCentralUnits() {
    near(central::lengthKm(1.0, 10.0), 14.76625, 1e-9, "1 M of a 10 solar-mass hole is 14.77 km");
    near(central::timeSeconds(1.0, 4.3e6), 21.18, 0.01, "1 M of time for Sgr A* is 21.2 s");
    // M87*: 1 M is 64 AU.
    near(central::lengthKm(1.0, 6.5e9) / central::kAuKm, 64.2, 0.1, "1 M of M87* is 64 AU");
    // Catalogue compactness R / M.
    for (const auto& e : central::catalog()) {
        if (std::string(e.label) == "Sun") near(e.object.radius, 4.711e5, 1e2, "the Sun has R = 4.71e5 M");
        if (std::string(e.label) == "Earth") near(e.object.radius / 1.44e9, 1.0, 0.01, "Earth has R = 1.44e9 M");
        if (std::string(e.label) == "Neutron star") near(e.object.radius, 5.8, 0.05, "a 12 km neutron star has R = 5.8 M");
    }
}

static void testSurfacePhysics() {
    // Neutron star, 1.4 solar masses, R = 5.8 M.
    near(1.0 + central::surfaceRedshift(5.8), 1.235441, 1e-6, "neutron star surface: 1 + z = 1.2355");
    near(central::surfaceGravity(5.8), 0.036726, 1e-6, "neutron star surface gravity 0.0367 / M");
    near(central::accelerationSI(central::surfaceGravity(5.8), 1.4) / 1.6e12, 1.0, 0.01,
         "neutron star surface gravity 1.6e12 m/s^2");
    // The Sun: z = 2.12e-6, measured as a 0.636 km/s redshift.
    const double rSun = central::kSunRadiusKm / central::kKmPerSolarMass;
    near(central::surfaceRedshift(rSun), 2.1225e-6, 2e-10, "solar surface redshift z = 2.12e-6");
    near(central::surfaceRedshift(rSun) * central::kSpeedOfLightKmS, 0.636, 0.001, "solar redshift is 0.636 km/s");
    // Sirius B at R = 3885 M: z = 2.575e-4 (77.2 km/s).
    near(central::surfaceRedshift(3885.0) * central::kSpeedOfLightKmS, 77.2, 0.1, "Sirius B redshift 77.2 km/s");
    // Limb deflection of starlight by the Sun: 4M/b = 1.751 arcsec, against the exact
    // integral (twice the sweep of a ray leaving the limb tangentially, minus pi).
    const double arcsec = kPi / (180.0 * 3600.0);
    const double b = rSun / std::sqrt(1.0 - 2.0 / rSun);
    near(4.0 / b / arcsec, 1.751, 0.001, "solar limb deflection 4M/b = 1.751 arcsec");
    near((2.0 * central::maxBendingAngle(rSun) - kPi) / (4.0 / b), 1.0, 1e-4, "exact solar deflection matches 4M/b");

    // A rigidly rotating surface element has a unit 4-velocity.
    {
        const double a = 0.3, R = 6.0, th = 1.0, omega = 0.02;
        const double ut = central::rotatingSurfaceUt(a, R, th, omega);
        // Kerr-Schild components equal Boyer-Lindquist ones when dr = 0.
        const Vec4 u = {ut, 0.0, 0.0, omega * ut};
        near(dot(metricTensor(a, R, th), u, u), -1.0, 1e-12, "rotating surface 4-velocity is unit");
    }
    near(central::Object{central::Kind::NeutronStar, "", 1.4, 5.8, 0.37, 0.35}.maxChi(), 0.548, 0.001,
         "mass-shedding limit of a neutron star: chi < 0.55");
}

// The angle a photon leaving radius R at angle alpha to the radial direction sweeps
// on its way to infinity, from the geodesic integrator (a = 0, equatorial).
static double tracedBending(double R, double alpha) {
    const double v[3] = {std::cos(alpha), 0.0, std::sin(alpha)};
    Constants c;
    State s;
    launchFromZamo(0.0, R, kPi / 2, 0.0, v, 0.0, c, s);
    while (s.r < 1e6) {
        stepMino(0.0, c, s, photonStepSize(0.0, c, s, 2e-4));
        enforceConstraints(0.0, c, s);
    }
    // What is left is a straight line in flat space: asin(b / r) more.
    return s.phi + std::asin(c.L / c.E / s.r);
}

static void testSurfaceLightBending() {
    // Neutron star, R = 5.8 M: exact integral against the ray tracer's integrator.
    for (double alphaDeg : {10.0, 40.0, 70.0, 89.0}) {
        const double alpha = alphaDeg * kPi / 180.0;
        near(tracedBending(5.8, alpha), central::bendingAngle(5.8, alpha), 1e-6,
             "traced light bending at R = 5.8 M matches the exact integral");
    }
    // Beloborodov's approximation, and how much of the star shows (~76%).
    const double psiMax = central::maxBendingAngle(5.8);
    near(psiMax * 180.0 / kPi, 121.75, 1.5, "light leaving the limb of R = 5.8 M sweeps ~122 degrees");
    near(central::bendingAngleBeloborodov(5.8, kPi / 2) * 180.0 / kPi, 121.75, 0.01, "Beloborodov: 121.75 degrees");
    near(central::visibleFraction(5.8), 0.76, 0.02, "a 5.8 M neutron star shows ~76% of its surface");
    const double rFull = central::fullVisibilityRadius();
    std::printf("      the whole surface is visible below R = %.4f M\n", rFull);
    near(rFull, 3.52, 0.02, "the whole surface is visible below R ~ 3.52 M");
    truth(central::visibleFraction(3.4) == 1.0, "R = 3.4 M: the whole surface shows");
}

static void testSurfaceRays() {
    // A radial photon from the static surface of a R = 5.8 M star reaches a static
    // camera far away with g = sqrt(1 - 2/R) / sqrt(1 - 2/r_cam).
    trajectory::RayScene scene;
    scene.disk = false;
    scene.surface = 5.8;
    scene.escapeRadius = 2e6;
    const double rCam = 1e6;
    const double atStar[3] = {-1.0, 0.0, 0.0};
    auto ray = trajectory::traceRay(0.0, rCam, kPi / 2, 0.0, zamoAt(0.0, rCam, kPi / 2), atStar, scene);
    truth(ray.end == trajectory::RayEnd::Surface, "looking at a neutron star: the ray ends on its surface");
    near(ray.g, std::sqrt(1.0 - 2.0 / 5.8) / std::sqrt(1.0 - 2.0 / rCam), 1e-9,
         "radial photon from a static surface: g = sqrt(1 - 2/R)");
    near(ray.emissionCos, 1.0, 1e-6, "... leaving along the normal");

    // The emission angle the ray tracer uses, cos = g |v_r| / sqrt(Sigma Delta), against
    // the surface element's own orthonormal frame (rotating Kerr surface).
    {
        const double a = 0.3, R = 6.0, th = 1.1, omega = 0.02;
        const double ut = central::rotatingSurfaceUt(a, R, th, omega);
        const Frame f = observerFrame(a, R, th, {ut, 0.0, 0.0, omega * ut});
        for (double angle : {0.2, 0.9, 1.4}) {
            // Leaving with unit energy along n = (cos, sin * 0.6, sin * 0.8) in that frame.
            const double n[3] = {std::cos(angle), 0.6 * std::sin(angle), 0.8 * std::sin(angle)};
            const double look[3] = {-n[0], -n[1], -n[2]};
            Constants c;
            State s;
            photonFromFrame(a, R, th, 0.3, f, look, c, s);
            const double g = 1.0 / (ut * (c.E - omega * c.L));
            const Metric m = metricAt(a, R, th);
            near(g, 1.0, 1e-12, "a photon of unit energy in the surface frame has g = 1");
            near(g * std::abs(s.vr) / std::sqrt(m.sigma * m.delta), n[0], 1e-12,
                 "emission angle from g |v_r| / sqrt(Sigma Delta) matches the surface frame");
        }
    }

    // Birkhoff: the same orbit around a black hole and a star; only how deep it reaches
    // differs. Periapsis 5 M, apoapsis 20 M (a = 0).
    double E = 0.0, L = 0.0;
    truth(equatorialOrbit(0.0, 5.0, 20.0, true, E, L), "bound orbit between 5 and 20 M exists");
    const Constants orbit = {E, L, 0.0, 1.0};
    const State start = {0.0, 20.0, kPi / 2, 0.0, 0.0, 0.0};
    auto around = [&](double surface) { return trajectory::classify(0.0, orbit, start, surface).fate; };
    truth(around(0.0) == trajectory::Fate::Bound, "around a black hole the orbit is bound");
    truth(around(4.5) == trajectory::Fate::Bound, "around a 4.5 M star the same orbit is bound (Birkhoff)");
    truth(around(5.8) == trajectory::Fate::HitsSurface, "around a 5.8 M star it hits the surface");
    const auto path = trajectory::predict(0.0, orbit, start, 2000.0, 0.3, 4000, 5.8);
    truth(path.fate == trajectory::Fate::HitsSurface && path.rMin >= 5.75, "predicted path ends on the surface");
}

static void testLens() {
    // A photon far from the hole (a = 0, r = 500) travelling along +phi passes a
    // small black hole (m = 0.01) 2 M closer in: it turns 4m/b = 0.02 rad toward it.
    const double r = 500.0, m = 0.01, b = 2.0;
    const double v[3] = {0.0, 0.0, 1.0};
    Constants c;
    State s;
    launchFromZamo(0.0, r, kPi / 2, 0.0, v, 0.0, c, s);
    const auto x = toCartesian(0.0, r, kPi / 2, 0.0);
    const trajectory::Point lens = {x[0] - b, x[1], x[2]};  // radially inward
    Constants forward = c;
    State sf = s;
    truth(trajectory::lensKick(0.0, lens, m, false, forward, sf), "photon passing at b = 2 M is not captured");
    const Metric mt = metricAt(0.0, r, kPi / 2);
    const double eLocal = forward.E / std::sqrt(mt.sigma * mt.delta / mt.A);
    near(sf.vr / (eLocal * std::sqrt(mt.sigma * mt.delta)), -std::sin(4.0 * m / b), 1e-12,
         "thin lens: the ray turns 4m/b toward the small black hole");
    near(eLocal, 1.0 / std::sqrt(1.0 - 2.0 / r), 1e-12, "... at fixed local frequency");
    near(radialPotential(0.0, forward, r) / (sf.vr * sf.vr), 1.0, 1e-9, "... and stays a null geodesic");
    Constants backward = c;
    State sb = s;
    trajectory::lensKick(0.0, lens, m, true, backward, sb);
    truth(sb.vr > 0.0, "traced backwards, the earlier direction is turned the other way");
    Constants captured = c;
    State sc = s;
    const trajectory::Point close = {x[0] - 0.05, x[1], x[2]};
    truth(!trajectory::lensKick(0.0, close, m, false, captured, sc), "b < 3 sqrt(3) m: captured");

    // Through the ray tracer: a sideways ray from r = 1000 passing a small black hole
    // (m = 0.05) at b = 3 leaves for the sky 4m/b = 0.067 rad off its unlensed course,
    // bent toward the lens.
    trajectory::RayScene scene;
    scene.disk = false;
    scene.escapeRadius = 5000.0;
    const double sideways[3] = {0.0, 0.0, 1.0};
    const Frame f = zamoAt(0.0, 1000.0, kPi / 2);
    auto plain = trajectory::traceRay(0.0, 1000.0, kPi / 2, 0.0, f, sideways, scene);
    scene.lenses.push_back({{1000.0, 400.0, 3.0}, 0.05});
    auto lensed = trajectory::traceRay(0.0, 1000.0, kPi / 2, 0.0, f, sideways, scene);
    truth(plain.end == trajectory::RayEnd::Sky && lensed.end == trajectory::RayEnd::Sky && lensed.lensings == 1,
          "the lensed ray still reaches the sky");
    const double turned = std::acos(std::clamp(plain.skyDirection[0] * lensed.skyDirection[0] +
                                               plain.skyDirection[1] * lensed.skyDirection[1] +
                                               plain.skyDirection[2] * lensed.skyDirection[2], -1.0, 1.0));
    near(turned, 4.0 * 0.05 / 3.0, 2e-3, "a small black hole deflects a passing ray by 4m/b");
    truth(lensed.skyDirection[2] > plain.skyDirection[2], "... toward itself");
}

static void testMercury() {
    // Mercury: a = 5.7909e10 m, e = 0.20563, 415.2 orbits per century: 42.98 arcsec.
    const double e = 0.20563;
    const double p = 5.7909e10 * (1.0 - e * e) / (central::kKmPerSolarMass * 1e3);
    const double perCentury = central::periapsisAdvance(p, e) * (36525.0 / 87.9691);
    near(perCentury * 180.0 / kPi * 3600.0, 42.98, 0.01, "Mercury's perihelion advances 42.98 arcsec per century");
}

static void testStellarStructure() {
    const central::TovStar star = central::solveTov(100.0, 2.0, 1.28e-3);
    std::printf("      TOV (K = 100, Gamma = 2, rho_c = 1.28e-3): M = %.4f, R = %.4f (%.2f km)\n", star.mass,
                star.radius, star.radius * central::kKmPerSolarMass);
    near(star.mass, 1.40, 0.005, "TOV polytrope: M = 1.40 solar masses");
    near(star.radius, 9.586, 0.01, "TOV polytrope: circumferential radius 9.586 (14.2 km)");
    const double mMax = central::neutronStarMaxMass();
    std::printf("      maximum mass %.4f\n", mMax);
    near(mMax, 1.637, 0.005, "TOV polytrope: maximum mass 1.64 solar masses");
    near(central::neutronStarRadiusKm(1.4), star.radius * central::kKmPerSolarMass, 0.05,
         "neutron-star radius for 1.4 solar masses from the mass-radius curve");
    near(central::whiteDwarfRadiusKm(0.6), 8657.0, 5.0, "Nauenberg: a 0.6 solar-mass white dwarf is 8,660 km");
    truth(central::whiteDwarfRadiusKm(1.5) == 0.0, "no white dwarf above the Chandrasekhar mass");
}

static void testAccretionEfficiency() {
    near(central::blackHoleEfficiency(0.0), 1.0 - std::sqrt(8.0 / 9.0), 1e-12, "black hole a = 0: 5.72%");
    near(central::blackHoleEfficiency(0.998), 0.321, 0.002, "black hole a = 0.998: 32%");
    const central::SurfaceAccretion ns = central::surfaceEfficiency(0.0, 5.8, 0.0);
    near(ns.rInner, 6.0, 1e-12, "neutron star R = 5.8: the disk stops at the ISCO");
    near(ns.total, 1.0 - std::sqrt(1.0 - 2.0 / 5.8), 1e-12, "neutron star R = 5.8: 19.1% in all");
    near(ns.disk, 1.0 - std::sqrt(8.0 / 9.0), 1e-12, "... of which the disk releases 5.72%");
    near(ns.boundaryLayer, std::sqrt(8.0 / 9.0) - std::sqrt(1.0 - 2.0 / 5.8), 1e-12, "... and the surface the rest");
    const central::SurfaceAccretion wd = central::surfaceEfficiency(0.0, 9800.0, 0.0);
    near(wd.total, 1.02e-4, 0.01e-4, "white dwarf: 1e-4 of the rest mass");
    near(central::alfvenRadiusKm(1e12, 12.0, 1.4, 1e17), 4225.0, 20.0, "Alfven radius for 1e12 G at 1e17 g/s");
}

static void testTidalScales() {
    near(central::hillsMass(central::kSunRadiusKm, 1.0) / 1e8, 1.14, 0.01, "Hills mass for a Sun-like star: 1.14e8");
    const double rStarSgr = central::kSunRadiusKm / (central::kKmPerSolarMass * 4.3e6);
    near(rStarSgr, 0.11, 0.001, "the Sun is 0.11 M across at Sgr A*");
    near(central::tidalRadius(rStarSgr, 1.0 / 4.3e6), 17.8, 0.1, "Sgr A* disrupts a Sun-like star at 18 M");
    const double rStarM87 = central::kSunRadiusKm / (central::kKmPerSolarMass * 6.5e9);
    truth(central::tidalRadius(rStarM87, 1.0 / 6.5e9) < horizonRadius(0.9), "M87* swallows a Sun-like star whole");
}

static void testCollapse() {
    std::string message;
    central::Object sun = central::catalog()[4].object;  // the Sun
    truth(std::string(sun.name) == "Sun", "catalogue entry 4 is the Sun");
    truth(!central::collapseIfUnstable(sun, message), "the Sun stands");
    sun.radius = 2.2;
    truth(central::collapseIfUnstable(sun, message) && sun.hasHorizon(), "below 9/4 M the Sun collapses");
    // The Sun's spin from its 25.4-day equatorial rotation and I = 0.07 M R^2.
    near(central::catalog()[4].object.chi, 0.22, 0.002, "the Sun's J/M^2 from its rotation period is 0.22");
    near(sun.chi, central::catalog()[4].object.chi, 1e-15, "... and it collapses keeping its spin");
    central::Object earth = central::catalog()[0].object;
    earth.radius = 2.0;
    truth(central::collapseIfUnstable(earth, message) && earth.chi == 0.998, "the Earth must shed spin to collapse");
    central::Object wd = central::catalog()[6].object;
    wd.massSolar = 1.5;
    truth(central::collapseIfUnstable(wd, message) && wd.kind == central::Kind::NeutronStar,
          "a white dwarf over 1.44 solar masses becomes a neutron star");
    central::Object ns = central::catalog()[8].object;
    ns.massSolar = 1.7;
    truth(central::collapseIfUnstable(ns, message) && ns.hasHorizon(), "a 1.7 solar-mass neutron star collapses");

    // The compactness slider's stages, and what the scene and the metric make of a body.
    truth(central::stageForCompactness(4.7e5, central::Kind::Star) == central::Kind::Star, "R = 4.7e5 M: a star");
    truth(central::stageForCompactness(1e4, central::Kind::Star) == central::Kind::WhiteDwarf, "R = 1e4 M: white dwarf");
    truth(central::stageForCompactness(10.0, central::Kind::Star) == central::Kind::NeutronStar, "R = 10 M: neutron star");
    central::Object squeezed = central::catalog()[4].object;
    squeezed.radius = 10.0;
    truth(central::becomeStage(squeezed, central::Kind::NeutronStar, 5772.0) &&
              squeezed.chi == central::catalog()[4].object.chi,
          "a squeezed Sun becomes neutron-star-like and keeps its J");
    near(central::catalog()[4].object.sceneScale(), 4.711e5 / 6.0, 50.0, "the Sun's scene is laid out in units of R/6");
    near(central::catalog()[8].object.sceneScale(), 1.0, 0.0, "a neutron star's scene is in M");
    near(central::catalog()[0].object.metricSpin(), 0.99, 0.0, "a planet's Kerr exterior uses a = 0.99, not 890");
    near(central::boundaryRadius(0.0, 0.0), 2.0, 1e-15, "a black hole's boundary is its horizon");
    near(central::displayTemperature(20000.0), 20000.0, 0.0, "below 30,000 K colors are drawn as they are");
    truth(central::displayTemperature(1e7) > central::displayTemperature(1e6) &&
              central::displayTemperature(1e7) < 2e5,
          "hotter neutron-star surfaces still look bluer, within the color table");
}

// --- Real objects at real sizes ------------------------------------------------------

// Periods between successive apoapses of the first body, from its r(t) sampled every
// `chunk` over `duration` (each maximum refined by a parabola through three samples).
static double apoapsisPeriod(BodySystem& sys, double spin, double chunk, double duration, double& rMin) {
    std::vector<double> times, radii;
    sys.setFrameTime(chunk, spin);
    for (double t = 0.0; t <= duration && !sys.bodies.empty(); t += chunk) {
        times.push_back(sys.time);
        radii.push_back(sys.bodies[0].state.r);
        sys.update(spin, chunk, false);
    }
    rMin = *std::min_element(radii.begin(), radii.end());
    std::vector<double> peaks;
    for (size_t i = 1; i + 1 < radii.size(); ++i) {
        if (radii[i] > radii[i - 1] && radii[i] >= radii[i + 1]) {
            const double a = radii[i - 1], b = radii[i], c = radii[i + 1];
            peaks.push_back(times[i] + 0.5 * chunk * (a - c) / (a - 2.0 * b + c));
        }
    }
    return peaks.size() >= 2 ? peaks[1] - peaks[0] : 0.0;
}

// A pulsar's light cylinder and beam (src/central.h): the geometry of the last open
// field line, against the measured relations.
static void testPulsarBeams() {
    constexpr double deg = 3.14159265358979 / 180.0;
    // c P / (2 pi): 1,608 km for the Crab's 33.7 ms.
    near(central::lightCylinderKm(0.0337), 299792.458 * 0.0337 / (2.0 * 3.14159265358979), 1e-9,
         "the Crab's light cylinder is c P / 2 pi (1,608 km)");
    near(central::radioEmissionHeightKm(1.0), 400.0, 1e-9, "radio comes from ~400 km above a 1 s pulsar");
    near(central::radioEmissionHeightKm(0.0337), 400.0 * std::pow(0.0337, 0.3), 1e-9,
         "and from 400 P^0.3 km (Kijak and Gil 2003): 145 km for the Crab");
    // Small angles: rho = 3/2 theta_open = 1.24 deg (r_em / 10 km)^(1/2) P^(-1/2) (Gil et al.).
    const double slow = 1.24 * deg * std::sqrt(central::radioEmissionHeightKm(1.0) / 10.0);
    near(central::beamHalfAngle(1.0) / slow, 1.0, 0.02, "a 1 s pulsar's beam: 1.24 deg sqrt(r_em / 10 km / P)");
    // Against Rankin's (1993) outer cones, 5.75 deg / sqrt(P): within the spread of
    // measured emission heights.
    near(central::beamHalfAngle(0.0337) / (5.75 * deg / std::sqrt(0.0337)), 1.0, 0.2,
         "the Crab's beam is ~28 deg wide, as the 5.75 deg / sqrt(P) of measured cones");
    truth(central::beamHalfAngle(7.5) < 6.0 * deg && central::beamHalfAngle(7.5) > 2.0 * deg,
          "a 7.5 s magnetar's beam is a few degrees");
    // Millisecond pulsars' open field lines fill most of the sky: very wide profiles.
    truth(central::beamHalfAngle(1.0 / 401.0) > 60.0 * deg, "a millisecond pulsar's beam is wider than 60 deg");
    truth(central::beamHalfAngle(0.1) > central::beamHalfAngle(1.0) && central::beamHalfAngle(1.0) > central::beamHalfAngle(10.0),
          "faster pulsars have wider beams");
    // A dipole line at colatitude theta points theta + atan(tan(theta) / 2) from the axis:
    // 90 deg where tan^2(theta) = 2; when more lines are open the beam fills the hemisphere.
    near(central::beamHalfAngle(1e-6), 0.5 * 3.14159265358979, 1e-12,
         "when every line is open the beam reaches the magnetic equator");

    // The core beam: 1.22 deg / sqrt(P) at half maximum (Rankin 1990), a fifth of the cone.
    near(central::beamCoreHalfWidth(1.0) / deg, 1.22, 1e-12, "a 1 s pulsar's core beam: 1.22 deg");
    near(central::beamCoreHalfWidth(0.0893) / deg, 4.08, 0.01, "Vela's core beam: 4.1 deg");
    // Torsional modes of the crust: f_l = f_2 sqrt((l - 1)(l + 2) / 4), f_2 = 22 Hz.
    near(central::crustTorsionalFrequency(2), 22.0, 1e-12, "the crust's fundamental torsional mode rings at 22 Hz");
    near(central::crustTorsionalFrequency(3) / central::crustTorsionalFrequency(2), std::sqrt(2.5), 1e-12,
         "l = 3 rings sqrt(5/2) times faster (35 Hz)");
    truth(central::crustTorsionalFrequency(2) > 18.0 && central::crustTorsionalFrequency(3) < 40.0,
          "within the 18-30+ Hz oscillations after SGR 1806-20's giant flare");

    // The catalogue's neutron star is a young pulsar like Vela: 89 ms, its beams drawn.
    const central::CatalogEntry* ns = nullptr;
    for (const central::CatalogEntry& e : central::catalog()) {
        if (std::string(e.label) == "Neutron star") ns = &e;
    }
    truth(ns != nullptr && ns->object.beacon && ns->object.pulsar && !ns->disk,
          "the neutron star is a young, isolated pulsar with its beams");
    if (ns) {
        const double period = 2.0 * 3.14159265358979 / ns->object.surfaceOmega() *
                              central::timeSeconds(1.0, ns->object.massSolar);
        near(period, 0.0893, 1e-6, "it turns in Vela's 89.3 ms");
        // Nearly aligned: its beams stay close to the spin axis while the star turns under
        // them, a full turn three times a second on screen.
        truth(ns->object.magneticTiltDeg <= 15.0, "its magnetic axis is nearly aligned with its spin");
        near(ns->timeScale * central::timeSeconds(1.0, ns->object.massSolar) / period, 3.0, 0.05,
             "three turns a second on screen");
    }
}

// A cataclysmic variable's geometry (src/central.h) against published fits.
static void testCloseBinary() {
    near(central::rocheLobeFraction(1.0), 0.3789, 0.0005, "Eggleton's Roche lobe for equal masses is 0.379 a");
    // L1 against Plavec and Kratochvil's (1964) fit, 0.500 - 0.227 log10(q), good to ~1%.
    bool l1 = true;
    for (double q : {0.1, 0.25, 0.5, 1.0}) {
        l1 = l1 && std::abs(central::innerLagrangeFraction(q) / (0.5 - 0.227 * std::log10(q)) - 1.0) < 0.015;
    }
    truth(l1, "L1 sits where Plavec and Kratochvil's fit puts it (q = 0.1 to 1)");
    near(central::innerLagrangeFraction(1.0), 0.5, 1e-9, "for equal masses L1 is halfway");
    // The stream's closest approach to the accretor, Lubow and Shu's (1975) fit
    // r_min / a = 0.0488 q^-0.464 (0.05 < q < 1).
    for (double q : {0.1, 0.25, 0.5}) {
        double closest = 0.0;
        central::rocheStream(q, 0.0, 16, closest);
        near(closest / (0.0488 * std::pow(q, -0.464)), 1.0, 0.06, "the stream swings within Lubow and Shu's r_min");
    }
    // Coriolis bends it ahead of the companion (toward +y) and it ends at the disk's rim.
    double closest = 0.0;
    const double rim = central::tidalTruncationFraction(0.25);
    const auto stream = central::rocheStream(0.25, rim, 16, closest);
    truth(stream.size() == 16 && stream.back()[1] > 0.0 && std::hypot(stream.back()[0], stream.back()[1]) < rim + 1e-3,
          "the stream bends ahead of the companion and stops at the disk's rim (0.48 a for q = 0.25)");
}

// The shape of a star filling its Roche lobe (the white dwarf's companion is drawn as it).
static void testRocheLobe() {
    // The potential's gradient (minus the effective gravity) against finite differences.
    const double q = 0.25, x = 0.7, y = 0.12, z = 0.05, h = 1e-6;
    const auto g = central::rocheGradient(q, x, y, z);
    auto phi = [&](double a, double b, double c) { return central::rochePotential(q, a, b, c); };
    const double fx = (phi(x + h, y, z) - phi(x - h, y, z)) / (2.0 * h);
    const double fy = (phi(x, y + h, z) - phi(x, y - h, z)) / (2.0 * h);
    const double fz = (phi(x, y, z + h) - phi(x, y, z - h)) / (2.0 * h);
    truth(std::abs(g[0] - fx) < 1e-6 && std::abs(g[1] - fy) < 1e-6 && std::abs(g[2] - fz) < 1e-6,
          "the Roche potential's gradient is its derivative");
    double force = 0.0, worst = 0.0, reachL1 = 0.0;
    bool teardrop = true;
    for (double qq : {0.1, 0.25, 1.0, 3.0}) {
        // No net force at L1 (a saddle of the potential).
        const double xL = central::innerLagrangeFraction(qq);
        force = std::max(force, std::abs(central::rocheGradient(qq, xL, 0.0, 0.0)[0]));
        // Eggleton (1983): the radius of a sphere of the lobe's volume is his fit's to 1%.
        const double rv = std::cbrt(3.0 * central::rocheLobeVolume(qq) / (4.0 * kerr::kPi));
        worst = std::max(worst, std::abs(rv / central::rocheLobeFraction(qq) - 1.0));
        // It reaches L1 toward the other star, and is a teardrop: shortest to the pole, wider
        // along the orbit, longer behind, longest toward L1 (Kopal 1959).
        const double point = central::rocheLobeReach(qq, -1.0, 0.0, 0.0);
        reachL1 = std::max(reachL1, std::abs(point - (1.0 - xL)));
        const double pole = central::rocheLobeReach(qq, 0.0, 0.0, 1.0);
        const double side = central::rocheLobeReach(qq, 0.0, 1.0, 0.0);
        const double back = central::rocheLobeReach(qq, 1.0, 0.0, 0.0);
        teardrop = teardrop && pole < side && side < back && back < point;
    }
    truth(force < 1e-9, "no net force at L1");
    near(worst, 0.0, 0.01, "Eggleton's radius matches the Roche lobe's volume to 1% (q = 0.1 to 3)");
    near(reachL1, 0.0, 1e-3, "a lobe-filling star reaches L1");
    truth(teardrop, "the lobe is a teardrop: pole < side < back < point (q = 0.1 to 3)");
}

static void testRealWorld() {
    using central::Host;
    const auto& cat = central::catalog();
    truth(cat.size() == 15, "the catalogue holds 15 objects");

    // Spins from observed rotation periods, chi = I Omega / M^2 (c J / G M^2).
    near(central::chiFromPeriod(0.0337, 1.4, 12.0, 0.35), 0.0152, 0.0002, "the Crab pulsar (33.7 ms): chi = 0.015");
    // J = 0.3307 M R^2 (2 pi / 86,164 s) = 5.86e33 kg m^2/s; G M^2 / c = 7.94e30 kg m^2/s.
    near(cat[0].object.chi, 5.86e33 / 7.94e30, 4.0, "the Earth's J/M^2 is ~738");
    near(cat[4].object.chi, 0.22, 0.002, "the Sun's (25.4 days) is 0.22");
    // The accreting millisecond pulsar's field holds its disk off at half the Alfven radius.
    near(0.5 * central::alfvenRadiusKm(2e8, 12.0, 1.4, 0.01 * central::eddingtonRate(1.4)), 26.0, 3.0,
         "a 2e8 G millisecond pulsar holds its disk off at ~26 km");
    truth(cat[9].object.pulsar && !cat[9].disk && cat[11].object.magnetar && !cat[11].disk,
          "the Crab pulsar and the magnetar have no disk");

    // Each catalogue object is the right kind of place.
    const Host want[] = {Host::Earthlike,   Host::GasGiant,    Host::BrownDwarf,  Host::RedDwarf,    Host::SunLike,
                         Host::Giant,       Host::WhiteDwarf,  Host::WhiteDwarf,  Host::NeutronStar, Host::NeutronStar,
                         Host::NeutronStar, Host::NeutronStar, Host::BlackHole,   Host::BlackHole,   Host::BlackHole};
    bool hosts = true;
    for (size_t k = 0; k < cat.size() && k < 15; ++k) hosts = hosts && central::hostOf(cat[k].object) == want[k];
    truth(hosts, "every catalogue object is classified as the right host");

    // Sizes in words.
    auto says = [](double km, const char* phrase) {
        return central::sizeComparison(km).find(phrase) != std::string::npos;
    };
    truth(says(12.0, "city"), "a neutron star is the size of a city");
    truth(says(8700.0, "1.4 times the size of the Earth"), "a white dwarf is 1.4 times the Earth");
    truth(says(central::kSunRadiusKm, "about the size of the Sun"), "the Sun is the size of the Sun");
    truth(says(100.0 * central::kSunRadiusKm, "Mercury's orbit"), "a red giant is as wide as Mercury's orbit");

    // Sunlight.
    near(central::fluxRelativeToEarth(central::kSunRadiusKm, central::kSunTemperature, central::kAuKm), 1.0, 1e-12,
         "the Sun's light at 1 AU is the unit");
    near(central::fluxRelativeToEarth(central::kSunRadiusKm, central::kSunTemperature, 0.5 * central::kAuKm), 4.0,
         1e-12, "... and falls off as 1/d^2");
    const central::Object& proxima = cat[3].object;
    const double rRatio = central::lengthKm(proxima.radius, proxima.massSolar) / central::kSunRadiusKm;
    near(rRatio * rRatio * std::pow(proxima.surfaceTemperature / central::kSunTemperature, 4.0) / 0.00165, 1.0, 0.15,
         "Proxima's radius and temperature give its luminosity (0.0017 Suns)");

    // Thin-disk temperatures: the peak factor at r = 49/36 r_in, by brute force.
    double best = 0.0;
    for (double x = 1.0; x < 4.0; x += 1e-5) best = std::max(best, std::pow(x, -0.75) * std::pow(1.0 - 1.0 / std::sqrt(x), 0.25));
    const double mdot = 0.01 * central::eddingtonRate(1.4), rIn = 12.4;
    const double tStar = std::pow(3.0 * 6.674e-8 * 1.4 * 1.98847e33 * mdot / (8.0 * kPi * 5.6704e-5 * std::pow(rIn * 1e5, 3.0)),
                                  0.25);
    near(central::diskPeakTemperature(1.4, mdot, rIn) / tStar, best, 1e-9, "the disk is hottest at r = 49/36 r_in");
    near(central::diskPeakTemperature(1.4, 2.0 * mdot, rIn) / central::diskPeakTemperature(1.4, mdot, rIn),
         std::pow(2.0, 0.25), 1e-12, "... and its temperature goes as mdot^(1/4)");
    const double tNs = central::diskPeakTemperature(1.4, mdot, rIn);
    const double tWd = central::diskPeakTemperature(0.6, 0.01 * central::eddingtonRate(0.6), 8700.0);
    truth(tNs > 1e6 && tNs < 1e7, "a neutron star's disk is X-ray hot (1e6-1e7 K)");
    truth(tWd > 1e4 && tWd < 5e4, "a white dwarf's disk is ultraviolet-hot (1e4-5e4 K)");

    // Bodies too small to resolve are drawn as points sending the same light.
    double dim = 0.0;
    const double drawn = pointSourceRadius(0.01, 100.0, 0.001, dim);
    near(drawn, kPointRadiusPixels * 0.1, 1e-15, "an unresolved body is drawn a pixel wide");
    near(drawn * drawn * dim, 1e-4, 1e-15, "... dimmed to send the same light");
    near(pointSourceRadius(1.0, 100.0, 0.001, dim), 1.0, 0.0, "a resolved body is drawn at its size");

    // What may orbit what.
    auto offered = [](Preset p, Host h) {
        for (const PresetInfo& info : presetList()) {
            if (info.preset == p) return presetAvailable(info, h);
        }
        return false;
    };
    truth(!offered(Preset::TidalDisruption, Host::NeutronStar) && !offered(Preset::Pulsar, Host::SunLike),
          "black-hole demonstrations are not offered around stars");
    truth(offered(Preset::Moon, Host::Earthlike) && offered(Preset::GalileanMoons, Host::GasGiant) &&
              offered(Preset::InnerPlanets, Host::SunLike) && offered(Preset::IscoOrbit, Host::NeutronStar),
          "real companions are offered where they belong");

    // Around a neutron star the orbit demonstrations carry a blob of plasma, 0.3 km.
    {
        BodySystem sys;
        const central::Object ns = cat[8].object;
        sys.host = central::hostOf(ns);
        sys.centralMassSolar = ns.massSolar;
        sys.surface = ns.radius;
        sys.setScale(ns.sceneScale(), ns.metricSpin());
        sys.addPreset(ns.metricSpin(), Preset::IscoOrbit);
        truth(sys.bodies.size() == 1 && sys.bodies[0].name == "Plasma blob", "a neutron star's ISCO holds plasma");
        if (!sys.bodies.empty()) {
            near(central::lengthKm(sys.bodies[0].radius, ns.massSolar), 0.3, 1e-6, "... 0.3 km across");
        }
    }

    // The Moon: real size, real orbit, and the period Kepler gives a test mass,
    // 2 pi sqrt(a^3 / G M_earth) = 27.45 days (the real 27.32 includes the Moon's own mass).
    {
        BodySystem sys;
        const central::Object earth = cat[0].object;
        const double spin = earth.metricSpin();
        sys.host = central::hostOf(earth);
        sys.centralMassSolar = earth.massSolar;
        sys.centralRadiusKm = central::kEarthRadiusKm;
        sys.surface = earth.radius;
        sys.setScale(earth.sceneScale(), spin);
        sys.addPreset(spin, Preset::Moon);
        truth(sys.bodies.size() == 1 && sys.bodies[0].kind == BodyKind::Moon, "the Moon is added");
        if (!sys.bodies.empty()) {
            const double M = earth.massSolar;
            near(central::lengthKm(sys.bodies[0].radius, M), 1737.4, 1e-3, "the Moon at its real radius");
            near(central::lengthKm(sys.bodies[0].state.r, M), 384400.0 * 1.0549, 1.0, "released at apogee, 405,500 km");
            const double a = central::kmToM(384400.0, M);
            const double kepler = 2.0 * kPi * a * std::sqrt(a);
            double rMin = 0.0;
            const double period = apoapsisPeriod(sys, spin, kepler / 2000.0, 2.2 * kepler, rMin);
            near(period / kepler, 1.0, 1e-3, "the Moon's orbit takes the Kepler period");
            near(central::timeSeconds(kepler, M) / 86400.0, 27.45, 0.05, "... 27.45 days");
            near(central::lengthKm(rMin, M), 384400.0 * (1.0 - 0.0549), 300.0, "perigee 363,300 km: the orbit holds");
        }
    }

    // Mercury round a Sun-like star, at 0.387 AU: its 88-day year.
    {
        BodySystem sys;
        const central::Object sun = cat[4].object;
        const double spin = sun.metricSpin();
        sys.host = central::hostOf(sun);
        sys.centralMassSolar = 1.0;
        sys.centralRadiusKm = central::kSunRadiusKm;
        sys.surface = sun.radius;
        sys.setScale(sun.sceneScale(), spin);
        sys.addPreset(spin, Preset::InnerPlanets);
        truth(sys.bodies.size() == 3 && sys.bodies[0].name == "Mercury", "Mercury, Venus and the Earth are added");
        if (sys.bodies.size() == 3) {
            near(central::lengthKm(sys.bodies[2].radius, 1.0), central::kEarthRadiusKm, 1e-3, "the Earth at its real size");
            const double a = central::kmToM(0.3871 * central::kAuKm, 1.0);
            const double kepler = 2.0 * kPi * a * std::sqrt(a);
            near(central::timeSeconds(kepler, 1.0) / 86400.0, 87.97, 0.1, "Mercury's year is 88 days");
            double rMin = 0.0;
            const double period = apoapsisPeriod(sys, spin, kepler / 2000.0, 2.2 * kepler, rMin);
            near(period / kepler, 1.0, 1e-3, "... and its orbit takes that long");
            near(central::lengthKm(rMin, 1.0) / central::kAuKm, 0.3871 * (1.0 - 0.2056), 0.001, "perihelion 0.307 AU");
        }
    }
}

// Runs one group of checks, and says so if it was slow.
// --- What surfaces look like ------------------------------------------------------

static void testSurfaceAppearance() {
    // Eddington-Barbier with a grey atmosphere: bolometrically, Eddington's law.
    for (double mu : {0.0, 0.3, 0.5, 1.0}) {
        const double ratio = std::pow(surface::limbTemperatureFactor(mu) / surface::limbTemperatureFactor(1.0), 4.0);
        near(ratio, (2.0 + 3.0 * mu) / 5.0, 1e-12, "grey limb darkening integrates to (2 + 3 mu) / 5");
    }
    // In visible light the limb is darker in the blue than in the red, and the Sun's
    // measured darkening at 550 nm (~0.70 at mu = 0.5, Neckel and Labs) is matched to ~5%.
    const double blue = surface::limbDarkening(5772.0, 0.2, 420.0), red = surface::limbDarkening(5772.0, 0.2, 680.0);
    truth(blue < red && red < 1.0, "the limb is darker and redder");
    near(surface::limbDarkening(5772.0, 0.5, 550.0), 0.70, 0.04, "solar limb darkening at 550 nm, mu = 0.5");

    // Convection: the Sun's granules (~1,000-1,500 km, 5-20 minutes); a red giant's
    // handful of huge cells.
    const double gSun = surface::surfaceGravityCgs(1.0, central::kSunRadiusKm);
    near(gSun, 27400.0, 100.0, "solar surface gravity (cm/s^2)");
    near(surface::pressureScaleHeightKm(5772.0, gSun), 134.0, 3.0, "solar pressure scale height (km)");
    const double granule = surface::granuleSizeKm(5772.0, gSun);
    truth(granule > 1000.0 && granule < 1500.0, "solar granules are ~1,300 km across");
    const double life = surface::granuleLifetimeSeconds(5772.0, gSun) / 60.0;
    truth(life > 5.0 && life < 20.0, "solar granules live 5-20 minutes");
    const double gGiant = surface::surfaceGravityCgs(1.0, 100.0 * central::kSunRadiusKm);
    const double cells = 100.0 * central::kSunRadiusKm / surface::granuleSizeKm(4000.0, gGiant);
    truth(cells > 3.0 && cells < 20.0, "a red giant has a few huge convection cells (per radian)");

    // Starspots: umbrae ~1,850 K cooler on the Sun, ~200 K on a cool red dwarf.
    near(surface::spotTemperatureDrop(5772.0), 1850.0, 1.0, "solar umbral temperature drop");
    const double proxima = surface::spotTemperatureDrop(3042.0);
    truth(proxima > 120.0 && proxima < 350.0, "red dwarf spots are only a few hundred K cooler");
    bool unit = true, belts = true;
    const surface::SpotPattern sun = surface::makeSpots(surface::SpotStyle::Solar, 0.6, 7u);
    for (const surface::Spot& spot : sun.spots) {
        unit = unit && std::abs(glm::length(spot.direction) - 1.0f) < 1e-5f;
        belts = belts && std::abs(std::asin(spot.direction.z)) < 45.0 * kPi / 180.0;
    }
    truth(!sun.spots.empty() && static_cast<int>(sun.spots.size()) <= surface::kMaxSpots, "spot count within the shader's list");
    truth(!sun.regions.empty() && static_cast<int>(sun.regions.size()) <= surface::kMaxRegions, "active regions");
    truth(unit, "spot directions are unit vectors");
    truth(belts, "sunspots stay in the low-latitude belts");
    const surface::SpotPattern dwarf = surface::makeSpots(surface::SpotStyle::RedDwarf, 0.8, 7u);
    double biggest = 0.0;
    for (const surface::Spot& spot : dwarf.spots) biggest = std::max(biggest, static_cast<double>(spot.radius));
    truth(biggest > 0.1, "red dwarfs carry far larger spots");

    // The corona: Allen's electron density times the dilution, integrated along the line
    // of sight, has the shape of Baumbach's measured brightness from 1.02 to 6 R.
    const double ref = surface::coronaColumn(1.1) / surface::baumbachBrightness(1.1);
    for (double rho : {1.02, 1.05, 1.2, 1.5, 2.0, 3.0, 6.0}) {
        const double shape = surface::coronaColumn(rho) / surface::baumbachBrightness(rho) / ref;
        near(shape, 1.0, 0.06, "K-corona profile against Baumbach");
    }
    near(surface::coronaScaleReal() * surface::coronaColumn(1.1), 0.8 * surface::baumbachBrightness(1.1), 1e-15,
         "corona calibrated to the disk center");

    // Which look each object gets.
    central::Object o = central::catalog()[4].object;  // the Sun
    surface::Appearance a = surface::appearance(o);
    truth(a.look == surface::Look::Photosphere && a.corona, "the Sun: granulation, spots and a corona");
    near(a.granuleFrequency, central::kSunRadiusKm / granule, 1e-6 * a.granuleFrequency, "granules per radian");
    truth(surface::appearance(central::catalog()[2].object).look == surface::Look::BrownDwarf, "a brown dwarf has weather");
    truth(surface::appearance(central::catalog()[0].object).atmosphere, "the Earth has air");
}

// --- Tidal disruption flares (src/tde.h) --------------------------------------------

// Coordinate time between two pericenter passages, by brute force: march in Mino time
// and time the moments vr turns from negative to positive. Also the first such moment.
static double measuredRadialPeriod(double a, const Constants& c, State s, double* firstPericenter = nullptr) {
    double prevVr = s.vr, prevT = s.t, first = -1.0;
    for (int i = 0; i < 40000000; ++i) {
        stepMino(a, c, s, 2e-5 * s.r);
        if (prevVr < 0.0 && s.vr >= 0.0) {
            const double f = -prevVr / (s.vr - prevVr);
            const double t = prevT + f * (s.t - prevT);
            if (first < 0.0) {
                first = t;
                if (firstPericenter) *firstPericenter = t;
            } else {
                return t - first;
            }
        }
        prevVr = s.vr;
        prevT = s.t;
    }
    return 0.0;
}

static void testTidalFlare() {
    using namespace tde;
    // A Sun torn apart by a million-Sun hole: the classic numbers. The most bound debris
    // (a = r_t^2 / 2 R*) returns after ~41 days; the fallback then peaks near 3 Suns a
    // year, ~100 times what the hole could swallow at its Eddington luminosity.
    const Scales sun6 = physicalScales(1e6, 1.0, central::kSunRadiusKm);
    near(sun6.tidalRadius / (central::kSunRadiusKm * 100.0), 1.0, 1e-12, "r_t = R* (M / m)^(1/3)");
    near(sun6.tMin / 86400.0, 40.95, 0.1, "first debris returns after 41 days (1e6 Suns)");
    near(sun6.peakRate, 2.97, 0.02, "peak fallback ~3 Suns a year");
    near(physicalScales(1e8, 1.0, central::kSunRadiusKm).tMin / sun6.tMin, 10.0, 1e-9, "t_min grows as M^(1/2)");
    const Scales g = classicalScales(1.0, 1e-3, 1.0);
    near(g.tMin, keplerPeriod(1.0, -g.energySpread), 1e-9 * g.tMin, "t_min is the period of the most bound debris");
    near(keplerPeriod(1.0, -1.0 / 200.0), 2.0 * kPi * std::pow(100.0, 1.5), 1e-9, "Kepler: 2 pi a^1.5");
    near(escapingLuminosity(0.5), 0.5, 1e-15, "below Eddington all the light escapes");
    near(escapingLuminosity(100.0), 1.0 + std::log(100.0), 1e-12, "photon trapping: 1 + ln(L / L_Edd)");
    near(eddingtonLuminosity(1e6) / 1.2572e44, 1.0, 1e-12, "Eddington luminosity of 1e6 Suns");

    // Radial periods from the quadrature against brute-force integration: Schwarzschild,
    // a rapidly spinning hole, and a very eccentric orbit like the debris.
    struct Case {
        double a, rp, ra;
        const char* what;
    };
    const Case cases[] = {{0.0, 8.0, 40.0, "radial period, Schwarzschild rp=8 ra=40"},
                          {0.9, 6.0, 30.0, "radial period, a=0.9 rp=6 ra=30"},
                          {0.6, 7.5, 120.0, "radial period, a=0.6 rp=7.5 ra=120 (eccentric)"}};
    for (const Case& k : cases) {
        double E, L;
        if (!equatorialOrbit(k.a, k.rp, k.ra, true, E, L)) {
            truth(false, k.what);
            continue;
        }
        const Constants c{E, L, 0.0, 1.0};
        // Start part way out, heading outward.
        const double r0 = 0.5 * (k.rp + k.ra);
        State s{0.0, r0, kPi / 2, 0.0, std::sqrt(radialPotential(k.a, c, r0)), 0.0};
        const RadialOrbit o = radialOrbit(k.a, c, s);
        truth(o.kind == RadialOrbit::Bound, "bound orbit recognised");
        near(o.rp, k.rp, 1e-7, "pericenter from the zeros of R(r)");
        near(o.ra, k.ra, 1e-6 * k.ra, "apocenter from the zeros of R(r)");
        double next = 0.0;
        const double measured = measuredRadialPeriod(k.a, c, s, &next);
        near(o.period / measured, 1.0, 2e-5, k.what);
        near(o.toPericenter / next, 1.0, 1e-4, "time to the next pericenter, outbound");
        State in = s;
        in.vr = -in.vr;
        const RadialOrbit back = radialOrbit(k.a, c, in);
        double nextIn = 0.0;
        measuredRadialPeriod(k.a, c, in, &nextIn);
        near(back.toPericenter / nextIn, 1.0, 1e-4, "time to the next pericenter, inbound");
    }
    {
        // Far out the period is Kepler's, with relativistic corrections of order M / a.
        double E, L;
        equatorialOrbit(0.0, 500.0, 1500.0, true, E, L);
        const Constants c{E, L, 0.0, 1.0};
        const State s{0.0, 500.0, kPi / 2, 0.0, 0.0, 0.0};
        near(radialOrbit(0.0, c, s).period / (2.0 * kPi * std::pow(1000.0, 1.5)), 1.0, 0.01, "wide orbit: Kepler period");
    }

    // A flat spread of orbital energies, -dE to +dE: half is bound, and the fallback
    // after the first return follows dM/dt = (1 / 2dE) (1/3) (2 pi)^(2/3) t^(-5/3).
    {
        const double dE = 0.03;
        std::vector<double> times, weights;
        const int n = 200000;
        for (int i = 0; i < n; ++i) {
            const double e = -dE + 2.0 * dE * (i + 0.5) / n;
            if (e < 0.0) {
                times.push_back(keplerPeriod(1.0, e));
                weights.push_back(1.0);
            }
        }
        Fallback f;
        f.build(times, weights, n);
        near(f.boundFraction, 0.5, 1e-6, "flat dM/dE: half the debris is bound");
        near(f.firstReturn / keplerPeriod(1.0, -dE), 1.0, 1e-4, "first return: the most bound debris");
        const double t = 3.0 * f.firstReturn;
        const double want = std::pow(2.0 * kPi, 2.0 / 3.0) / (3.0 * 2.0 * dE) * std::pow(t, -5.0 / 3.0);
        near(f.rate(t) / want, 1.0, 0.02, "fallback rate at 3 t_min against the t^(-5/3) law");
        near(f.lateSlope, -5.0 / 3.0, 0.03, "fallback fades as t^(-5/3)");
        near(f.returned(1e12), 0.5, 0.01, "all the bound half comes back");
        truth(f.peakTime > f.firstReturn && f.peakTime < 1.6 * f.firstReturn, "peak shortly after the first return");
    }

    // The wind: closed-form travel time against a direct integral of dR / v(R).
    {
        const double vInf = 0.1, R0 = 15.0, R1 = 300.0;
        const int n = 20000;
        double sum = 0.0;
        for (int i = 0; i <= n; ++i) {
            const double R = R0 + (R1 - R0) * i / n;
            sum += (i == 0 || i == n ? 1.0 : (i % 2 ? 4.0 : 2.0)) / windSpeed(1.0, vInf, R);
        }
        sum *= (R1 - R0) / n / 3.0;
        near(windTravelTime(1.0, vInf, R0, R1) / sum, 1.0, 1e-8, "wind travel time in closed form");
        near(windTravelTime(1.0, vInf, R0, windFront(1.0, vInf, R0, 1234.5)), 1234.5, 1e-6, "wind front inverts it");
        near(windSpeed(1.0, vInf, 1e12), vInf, 1e-6, "the wind coasts at its terminal speed");
        near(openingAngle(0.3, 0.9, 0.0), 0.9, 1e-15, "wide outflow before a disk forms");
        near(openingAngle(0.3, 0.9, 0.2), 0.3 + 0.6 / std::exp(1.0), 1e-12, "collimation e-folds at a fifth returned");
        near(openingAngle(0.3, 0.9, 5.0), 0.3, 1e-9, "narrow funnel once the disk has formed");
    }

    // A real disruption: a star 1.5 M in radius torn apart at pericenter (7.5 M) of a
    // parabolic orbit around a hole of spin 0.7. The energy spread straddles E = 1:
    // about half escapes, and the rest returns after roughly the period of the most
    // bound debris, with a spread dE ~ M R* / r_p^2.
    {
        const double a = 0.7, rp = 7.5, rStar = 1.5;
        const double A = rp - 2.0, B = 4.0 * a, C = -2.0 * (a * a + rp * rp);
        const double L = (-B + std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
        const Constants c{1.0, L, 0.0, 1.0};
        truth(std::abs(radialPotential(a, c, rp)) < 1e-9, "parabolic orbit turns at its pericenter");
        double v[3];
        zamoVelocity(a, c, State{0.0, rp, kPi / 2, 0.0, 0.0, 0.0}, v);
        std::mt19937 rng(11);
        const tidal::Mat3 F = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
        const tidal::Mat3 Fdot = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
        const auto gas = debris::deformedStar(a, rp, kPi / 2, 0.0, v, rStar, F, Fdot, 40000, 9000.0f, rng);
        const auto start = std::chrono::steady_clock::now();
        const Event ev = analyze(a, 100.0, gas);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        truth(seconds < 1.0, "analysing a disruption takes well under a second");
        near(ev.unboundFraction, 0.5, 0.08, "parabolic disruption: about half the debris escapes");
        const double tMin = keplerPeriod(1.0, -rStar / (rp * rp));
        truth(ev.fallback.firstReturn > 0.3 * tMin && ev.fallback.firstReturn < 1.5 * tMin,
              "first return near the period of the most bound debris");
        truth(ev.fallback.peakTime > ev.fallback.firstReturn, "fallback peaks after the first return");
        near(ev.fallback.lateSlope, -5.0 / 3.0, 0.35, "fallback of real debris fades roughly as t^(-5/3)");
    }

    // Labels across the stretch: a star stretched along x (F = diag(5, 1, 1/2)) is
    // labelled by where each parcel sat in y and z.
    {
        std::mt19937 rng(5);
        const tidal::Mat3 F = {{{5.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}}};
        const tidal::Mat3 Fdot = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
        const double still[3] = {0.0, 0.0, 0.0};
        const auto gas = debris::deformedStar(0.0, 50.0, kPi / 2, 0.0, still, 1.0, F, Fdot, 3000, 5000.0f, rng);
        double worst = 0.0;
        for (const ParticleInit& p : gas) {
            const double y = p.r * std::sin(p.theta) * std::sin(p.phi), z = p.r * std::cos(p.theta);
            const double want = y * y + (z / 0.5) * (z / 0.5);
            worst = std::max(worst, std::abs(p.label[0] * p.label[0] + p.label[1] * p.label[1] - want));
        }
        near(worst, 0.0, 2e-4, "debris labels: position across the stretch, in star radii");
    }
}

static void testMerger() {
    using namespace merger;
    using merger::kPi;
    // GW170817: two stars of ~1.365 Suns have a chirp mass of 1.188 Suns, and their chirp
    // took ~100 s to sweep up from 24 Hz to the merger.
    near(chirpMass(1.365, 1.365), 1.1883, 1e-3, "GW170817 chirp mass 1.188 Suns");
    near(chirpFrequency(101.0 / kSecond, chirpMass(1.365, 1.365)) / kSecond, 24.0, 0.5,
         "GW170817 was at 24 Hz ~100 s before the merger");
    {
        // The Peters time and the chirp agree: f = 2 f_orbit = sqrt(M / a^3) / pi.
        const double a = 40.0, m = 1.35;
        near(chirpFrequency(coalescenceTime(a, m, m), chirpMass(m, m)) / (std::sqrt(2 * m / (a * a * a)) / kPi), 1.0,
             1e-9, "chirp frequency = twice the Kepler frequency at the Peters separation");
        near(separationAfter(a, 0.5 * coalescenceTime(a, m, m), m, m), a * std::pow(0.5, 0.25), 1e-9,
             "a(t) = a0 (1 - t/T)^(1/4)");
    }
    {
        // The 2.5PN reaction force drains a circular orbit at Peters' rate, and integrating
        // it (with Newtonian gravity) reproduces his inspiral time.
        const double m1 = 1.35, m2 = 1.25, M = m1 + m2, mu = m1 * m2 / M;
        const double a0 = 40.0;
        glm::dvec3 x(a0, 0, 0), v(0, std::sqrt(M / a0), 0);
        near(mu * glm::dot(v, radiationReaction(x, v, m1, m2)) / -gwLuminosity(a0, m1, m2), 1.0, 1e-12,
             "radiation reaction loses energy at the quadrupole rate");
        auto accel = [&](const glm::dvec3& p, const glm::dvec3& w) {
            const double r = glm::length(p);
            return -M * p / (r * r * r) + radiationReaction(p, w, m1, m2);
        };
        double t = 0.0;
        while (glm::length(x) > 12.0) {
            const double r = glm::length(x);
            const double h = 0.002 * std::sqrt(r * r * r / M);
            const glm::dvec3 k1x = v, k1v = accel(x, v);
            const glm::dvec3 k2x = v + 0.5 * h * k1v, k2v = accel(x + 0.5 * h * k1x, v + 0.5 * h * k1v);
            const glm::dvec3 k3x = v + 0.5 * h * k2v, k3v = accel(x + 0.5 * h * k2x, v + 0.5 * h * k2v);
            const glm::dvec3 k4x = v + h * k3v, k4v = accel(x + h * k3x, v + h * k3v);
            x += h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
            v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
            t += h;
        }
        const double peters = coalescenceTime(a0, m1, m2) - coalescenceTime(12.0, m1, m2);
        near(t / peters, 1.0, 0.01, "integrated inspiral from 40 M to 12 M takes Peters' time");
    }
    {
        // Quadrupole waves of a circular binary: amplitude 4 mu M / (a D) face-on, half that
        // in h+ (and no hx) edge-on.
        const double m1 = 1.4, m2 = 1.2, M = m1 + m2, mu = m1 * m2 / M, a = 30.0, D = 1e6;
        const double w = std::sqrt(M / (a * a * a));
        for (double phase : {0.0, 0.7, 2.1}) {
            const glm::dvec3 n(std::cos(phase), std::sin(phase), 0.0);
            const glm::dvec3 t(-n.y, n.x, 0.0);
            std::vector<glm::dvec3> xs = {m2 / M * a * n, -m1 / M * a * n};
            std::vector<glm::dvec3> vs = {m2 / M * a * w * t, -m1 / M * a * w * t};
            std::vector<glm::dvec3> as = {-m2 / M * a * w * w * n, m1 / M * a * w * w * n};
            const auto I = quadrupoleSecondDerivative(xs, vs, as, {m1, m2});
            const Strain faceOn = strain(I, glm::dvec3(0, 0, 1), D);
            near(std::hypot(faceOn.plus, faceOn.cross) / (4.0 * mu * M / (a * D)), 1.0, 1e-9,
                 "face-on strain 4 mu M / (a D)");
            const Strain edgeOn = strain(I, glm::dvec3(0, 1, 0), D);
            near(edgeOn.cross * D, 0.0, 1e-9, "edge-on: no cross polarization");
            near(std::abs(edgeOn.plus) / (2.0 * mu * M / (a * D) * std::abs(std::cos(2.0 * phase))), 1.0, 1e-9,
                 "edge-on: h+ at half the face-on amplitude");
        }
    }
    {
        // An n = 1 polytrope of 1.35 Suns and 12 km: central density 1.2e15 g/cm^3, in
        // hydrostatic equilibrium, dP/dr = -rho m / r^2.
        const double M = 1.35, R = 12.0 / kKm;
        near(polytropeCentralDensity(M, R) * kGramsPerCc / 1.22e15, 1.0, 0.01, "central density 1.2e15 g/cc");
        near(polytropeEnclosed(M, R, R), M, 1e-12, "the polytrope's mass");
        const double K = polytropeK(R), r = 0.5 * R, dr = 1e-5;
        const double dP = K * (std::pow(polytropeDensity(M, R, r + dr), 2) - std::pow(polytropeDensity(M, R, r - dr), 2)) /
                          (2 * dr);
        near(dP / (-polytropeDensity(M, R, r) * polytropeEnclosed(M, R, r) / (r * r)), 1.0, 1e-6,
             "polytrope in hydrostatic equilibrium");
        const auto ps = starParticles(5000, R);
        near(static_cast<double>(ps.size()), 5000.0, 0.0, "star built from exactly n particles");
        double rMax = 0.0;
        int inner = 0;
        for (const auto& p : ps) {
            rMax = std::max(rMax, glm::length(p));
            if (glm::length(p) < 0.5 * R) ++inner;
        }
        truth(rMax < R, "particles inside the star");
        near(inner / 5000.0, polytropeEnclosed(1.0, R, 0.5 * R), 2e-3, "particles follow the polytrope's mass profile");
        near(temperatureFromThermalEnergy(thermalEnergyFromTemperature(3e10)), 3e10, 1.0, "temperature round trip");
        // 10 MeV ~ 1.16e11 K is u_th ~ 0.016 for an ideal nucleon gas.
        near(thermalEnergyFromTemperature(1.1605e11), 1.5 * 10.0 / 939.565, 1e-4, "kT = 10 MeV");
    }
    {
        // Collapse: a cold ball falls in whole; a Keplerian ring far out stays behind.
        std::vector<glm::dvec3> x, v;
        std::vector<double> m, phi;
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        while (x.size() < 2000) {
            glm::dvec3 p(u(rng), u(rng), u(rng));
            if (glm::dot(p, p) > 1.0) continue;
            x.push_back(4.0 * p);
            v.push_back(glm::dvec3(0.0));
            m.push_back(2.5 / 2000.0);
            phi.push_back(-1.0);
        }
        for (int k = 0; k < 200; ++k) {
            const double ang = 2.0 * kPi * k / 200.0, r = 40.0, vc = std::sqrt(2.5 / r);
            x.push_back({r * std::cos(ang), r * std::sin(ang), 0.0});
            v.push_back({-vc * std::sin(ang), vc * std::cos(ang), 0.0});
            m.push_back(0.1 / 200.0);
            phi.push_back(-2.5 / r);
        }
        std::vector<bool> alive(x.size(), true), swallowed;
        const BlackHole bh = collapse(x, v, m, phi, alive, glm::dvec3(0.0), swallowed);
        near(bh.mass, 2.5, 1e-9, "the cold ball collapses whole");
        near(bh.diskMass, 0.1, 1e-9, "the ring at 40 M stays out as a disk");
        near(bh.spin, 0.0, 1e-9, "a non-rotating ball makes a non-rotating hole");
    }
    {
        // Ringdown (Berti et al. 2006): M omega = 0.3737 for Schwarzschild; a 2.6-Sun hole
        // of spin 0.7 rings at ~6.6 kHz.
        near(ringdown(1.0, 0.0).omega, 0.3737, 0.01, "Schwarzschild l=m=2 QNM frequency");
        near(ringdown(1.0, 0.0).tau, 11.24, 0.6, "Schwarzschild QNM damping time ~11 M");
        near(ringdown(2.6, 0.7).omega / (2 * kPi) / kSecond / 1000.0, 6.65, 0.1, "2.6 Suns, chi 0.7: rings at 6.6 kHz");
        near(horizonRadius(2.0, 0.0), 4.0, 1e-12, "Schwarzschild horizon 2M");
    }
    {
        // r-process heating ~2e10 erg/g/s a day after the merger (Korobkin et al. 2012;
        // Metzger et al. 2010), falling as t^-1.3.
        near(rProcessHeating(kDay) / 2.0e10, 1.0, 0.1, "r-process heating 2e10 erg/g/s at 1 day");
        near(std::log(rProcessHeating(10 * kDay) / rProcessHeating(kDay)) / std::log(10.0), -1.3, 0.01,
             "heating falls as t^-1.3");
        const auto ej = gw170817Ejecta();
        const Kilonova blue(ej[0]), red(ej[2]);
        truth(blue.peakTime() > 0.3 * kDay && blue.peakTime() < 2.0 * kDay, "blue kilonova peaks within ~a day");
        truth(red.peakTime() > blue.peakTime(), "the red component peaks later");
        double total = 0.0;
        for (const auto& c : ej) total += Kilonova(c).luminosity(kDay);
        truth(total > 2e41 && total < 3e42, "AT2017gfo: ~5e41 erg/s a day after the merger");
        truth(blue.temperature(0.5 * kDay) > 6000.0, "blue component hot at half a day");
        near(red.temperature(10 * kDay), 3745.0, 1.0, "red photosphere held at its floor late on");
    }
    {
        // The debris cloud the collision throws out: at 0.3 c its edge is 90 km further out
        // each millisecond; it cools adiabatically as 1 / R (radiation dominated) and its
        // column thins as 1 / R^2; a few thousandths of a Sun are opaque for hours even to
        // electron scattering alone.
        const DebrisCloud cloud;
        near(debrisRadiusKm(cloud, 0.0), cloud.launchKm, 1e-12, "the cloud starts at the stars' surface");
        near(debrisRadiusKm(cloud, 1e-3) - debrisRadiusKm(cloud, 0.0), 89.94, 0.01, "0.3 c: 90 km per ms");
        near(debrisTemperature(cloud, 1e-3) * debrisRadiusKm(cloud, 1e-3), cloud.launchTemperature * cloud.launchKm,
             1e-3 * cloud.launchTemperature * cloud.launchKm, "T R stays constant as it expands");
        near(debrisOpticalDepth(cloud, 2e-3) / debrisOpticalDepth(cloud, 1e-2),
             std::pow(debrisRadiusKm(cloud, 1e-2) / debrisRadiusKm(cloud, 2e-3), 2.0), 1e-9,
             "its optical depth falls as 1 / R^2");
        truth(debrisOpticalDepth(cloud, 1e-3) > 1e12, "a millisecond after the merger it is utterly opaque");
        truth(debrisOpticalDepth(cloud, 3600.0) > 100.0, "still opaque an hour later");
    }
    {
        // GRB 170817A: gamma-rays 1.7 s after the merger, from a jet that broke out of the
        // ejecta; a jet core seen 25 degrees off-axis is de-boosted ~10 times below 1.
        near(breakoutTime(Jet{}, kEjectaEdge * kBluePolarStretch * gw170817Ejecta()[0].velocity), 1.7, 0.02,
             "jet breaks out of the ejecta 1.7 s after the merger");
        near(dopplerFactor(100.0, 0.0), 2.0 * 100.0, 0.02 * 200.0, "on-axis Doppler factor ~2 Gamma");
        truth(dopplerFactor(100.0, 25.0 * kPi / 180.0) < 0.15, "seen off-axis the core is de-boosted");
        const Tempo tp;
        near(tempo(tp, -1.0, 0.02, 1.7), tp.inspiralRate, 1e-12, "one second before: 0.32 s per s");
        near(tempo(tp, -1e-7, 0.02, 1.7), tp.mergerRate, 1e-12, "at the merger: slow motion");
    }
}

// --- Saving and loading scenes (src/scene_file.h, src/json.h) ---------------------------

// A scene as the scene file writes it.
static std::string sceneText(const Settings& settings, const Camera& camera, const BodySystem& bodies,
                             std::vector<char>* binary = nullptr) {
    scene::Blob blob;
    const std::string text = json::dump(scene::save(settings, camera, bodies, blob));
    if (binary) *binary = blob.serialize();
    return text;
}

static void testSceneFile() {
    {
        // Numbers come back bit for bit, the hard ones included (the smallest denormal,
        // the largest double, a negative zero; floats in their own shortest form).
        const double doubles[] = {0.1, 1.0 / 3.0, 5e-324, 1.7976931348623157e308, -2.5e-310, 6.02214076e23, -0.0,
                                  kPi};
        const float floats[] = {0.1f, 0.05f, 1e-45f, 3.4028235e38f, 15.0f, 6500.0f, 0.94f};
        json::Value list = json::Value::array();
        for (double d : doubles) list.push(json::Value::number(d));
        for (float f : floats) list.push(json::Value::number(f));
        list.push(json::Value::number(INFINITY));
        list.push(json::Value::number(-INFINITY));
        list.push(json::Value::number(NAN));
        const std::string odd = "quote \" backslash \\ line\nbreak, tab\t, \xC3\xA9, \xE2\x9C\x93, \x01";
        list.push(json::Value::string(odd));
        json::Value back;
        std::string error;
        const bool parsed = json::parse(json::dump(list), back, error);
        truth(parsed, "JSON: what is written parses");
        bool exact = parsed && back.size() == list.size();
        size_t i = 0;
        for (double d : doubles) {
            double x = 1.0;
            exact = exact && back[i++].get(x) && std::memcmp(&x, &d, sizeof d) == 0;
        }
        for (float f : floats) {
            float x = 1.0f;
            exact = exact && back[i++].get(x) && std::memcmp(&x, &f, sizeof f) == 0;
        }
        truth(exact, "JSON: doubles and floats come back bit for bit");
        double inf = 0.0, minusInf = 0.0, nan = 0.0;
        truth(parsed && back[i].get(inf) && back[i + 1].get(minusInf) && back[i + 2].get(nan) && inf == INFINITY &&
                  minusInf == -INFINITY && std::isnan(nan),
              "JSON: infinities and NaN survive (as strings)");
        std::string text;
        truth(parsed && back[i + 3].get(text) && text == odd, "JSON: escapes and UTF-8 survive");
        json::Value bad;
        const bool refused = !json::parse("{\n  \"a\": [1, 2,\n}", bad, error);
        truth(refused && error.find("line 3") != std::string::npos, "JSON: a malformed file is refused, with its line");
        truth(!json::parse("[1e999]", bad, error), "JSON: an out-of-range number is refused");
    }

    // A scene with everything in it: a star on its way to being torn apart, a binary, a
    // precessing orbit, bodies pulling on each other, a camera falling in.
    const double spin = 0.7;
    Settings settings;
    settings.central.chi = spin;
    settings.spin = static_cast<float>(spin);
    settings.diskOuter = 17.3f;
    settings.timeScale = 33.3f;
    settings.debrisFilaments = true;
    Camera camera;
    camera.motion = CameraMotion::FreeFall;
    {
        const double v[3] = {-0.1, 0.02, 0.05};
        launchFromZamo(spin, 30.0, 1.2, 0.4, v, 1.0, camera.geodesic, camera.path);
        camera.r = camera.path.r;
        camera.theta = camera.path.theta;
        camera.phi = camera.path.phi;
        camera.properTime = 3.25;
    }
    BodySystem original;
    original.addPreset(spin, Preset::TidalDisruption);
    original.addPreset(spin, Preset::BinaryStar);
    original.addPreset(spin, Preset::Precessing);
    auto runFor = [&](BodySystem& b, int steps) {
        for (int k = 0; k < steps; ++k) {
            b.setFrameTime(1.0, spin);
            b.update(spin, 1.0, true);
            b.pendingDebris.clear();  // the app hands these to the GPU every frame
        }
    };
    runFor(original, 450);
    truth(!original.tidalEvents.empty() && !original.departed.empty(), "test scene: the star has been torn apart");

    std::vector<char> binary1;
    const std::string text1 = sceneText(settings, camera, original, &binary1);
    json::Value root;
    scene::Blob blob;
    std::string error;
    const bool read = json::parse(text1, root, error) && blob.deserialize(binary1, error);
    Settings settings2;
    Camera camera2;
    BodySystem copy;
    const bool loaded = read && scene::load(root, blob, settings2, camera2, copy, error);
    truth(loaded, "scene: loads what was saved");
    std::vector<char> binary2;
    const std::string text2 = sceneText(settings2, camera2, copy, &binary2);
    truth(text1 == text2 && binary1 == binary2, "scene: save, load, save gives identical files");

    // ... and it runs on exactly as the original does, random numbers included.
    runFor(original, 150);
    runFor(copy, 150);
    original.addPreset(spin, Preset::Cluster);
    copy.addPreset(spin, Preset::Cluster);
    runFor(original, 50);
    runFor(copy, 50);
    std::vector<char> after1, after2;
    truth(sceneText(settings, camera, original, &after1) == sceneText(settings2, camera2, copy, &after2) &&
              after1 == after2,
          "scene: the loaded copy evolves bit for bit like the original");

    // Through files on disk; a sidecar that does not match is caught.
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "bh_test_scene.bhscene";
    const bool wrote = scene::writeFiles(path, root, blob, error);
    json::Value root3;
    scene::Blob blob3;
    const bool reread = wrote && scene::readFiles(path, root3, blob3, error);
    truth(reread && blob3.serialize() == binary1, "scene: written to disk and read back");
    {
        std::fstream f(scene::sidecarPath(path), std::ios::in | std::ios::out | std::ios::binary);
        f.seekp(100);
        f.put('\x7f');
    }
    truth(!scene::readFiles(path, root3, blob3, error) && error.find("does not belong") != std::string::npos,
          "scene: a changed sidecar is refused");
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(scene::sidecarPath(path), ec);

    // A bad value is refused, says where it is, and changes nothing.
    {
        std::string broken = text1;
        const size_t at = broken.find("\"kind\": 0");
        broken.replace(at, 9, "\"kind\": 99");
        json::Value badRoot;
        Settings untouched;
        untouched.diskOuter = 5.0f;
        Camera cam;
        BodySystem bodies;
        const bool parsedBroken = json::parse(broken, badRoot, error);
        const bool refused = parsedBroken && !scene::load(badRoot, blob, untouched, cam, bodies, error);
        truth(refused && error.find("settings.central.kind") != std::string::npos && untouched.diskOuter == 5.0f &&
                  bodies.bodies.empty(),
              "scene: an impossible value is refused, with where it is, and nothing changes");
        std::string newer = text1;
        newer.replace(newer.find("\"version\": 1"), 12, "\"version\": 99");
        truth(json::parse(newer, badRoot, error) && !scene::load(badRoot, blob, untouched, cam, bodies, error) &&
                  error.find("newer") != std::string::npos,
              "scene: a file from a newer version is refused");
    }
}

// --- The timeline (src/timeline.h) ------------------------------------------------------

static void testTimeline() {
    // A scene of bodies pulling on each other, stepped by frames of uneven length (as
    // real frames are), with checkpoints every 25 M and the user adding bodies at step 150.
    const double spin = 0.6;
    const Settings settings{};
    const Camera camera{};
    auto frameDt = [](size_t i) { return 0.6 + 0.35 * std::sin(0.37 * static_cast<double>(i)); };
    auto stepScene = [&](BodySystem& b, double dt) {
        b.setFrameTime(dt, spin);
        b.update(spin, dt, true);
        b.pendingDebris.clear();
    };
    auto checkpointOf = [&](const BodySystem& b, bool pinned) {
        Timeline::Checkpoint c;
        c.pinned = pinned;
        c.scene = sceneText(settings, camera, b);
        scene::Blob blob;
        scene::save(settings, camera, b, blob);
        c.blob = blob;
        return c;
    };
    auto restore = [&](const Timeline::Checkpoint& c, BodySystem& b) {
        json::Value root;
        std::string error;
        Settings s;
        Camera cam;
        return json::parse(c.scene, root, error) && scene::load(root, c.blob, s, cam, b, error);
    };
    const size_t kSteps = 400, kAction = 150;
    BodySystem live;
    live.addPreset(spin, Preset::BinaryStar);
    live.addPreset(spin, Preset::PlanetWithMoon);
    Timeline timeline;
    timeline.addCheckpoint(checkpointOf(live, true), live.time);
    std::vector<std::string> expected(kSteps + 1);  // the scene after every step, for comparison
    expected[0] = sceneText(settings, camera, live);
    for (size_t i = 0; i < kSteps; ++i) {
        if (i == kAction) {
            live.addPreset(spin, Preset::Precessing);
            timeline.addCheckpoint(checkpointOf(live, true), live.time);
            expected[i] = sceneText(settings, camera, live);
        }
        const double dt = frameDt(i);
        stepScene(live, dt);
        timeline.addStep({dt, 1.0f / 60.0f, {}}, live.time);
        expected[i + 1] = sceneText(settings, camera, live);
        if (live.time - timeline.lastCheckpointTime() >= 25.0) timeline.addCheckpoint(checkpointOf(live, false), live.time);
    }

    // Seeking: restore the checkpoint before, replay the logged steps.
    auto seek = [&](size_t target, BodySystem& b) {
        const Timeline::Checkpoint& c = timeline.checkpointFor(target);
        if (!restore(c, b)) return false;
        for (size_t i = c.step; i < target; ++i) stepScene(b, timeline.step(i).dt);
        timeline.setPosition(target);
        return true;
    };
    bool exact = true;
    for (size_t target : {size_t(0), size_t(37), size_t(149), size_t(150), size_t(151), size_t(288), kSteps}) {
        BodySystem b;
        exact = exact && seek(target, b) && sceneText(settings, camera, b) == expected[target];
    }
    truth(exact, "timeline: any moment, restored and replayed, is the scene as it was, bit for bit");
    truth(timeline.checkpointFor(kAction + 1).step >= kAction && timeline.checkpointFor(kAction).step == kAction,
          "timeline: no replay runs across what the user did");
    const size_t at = timeline.positionAt(137.3);
    truth(timeline.timeAt(at) <= 137.3 && timeline.timeAt(at + 1) > 137.3, "timeline: finds the frame at a time");
    truth(timeline.positionAt(-5.0) == 0 && timeline.positionAt(1e9) == kSteps, "timeline: times outside it clamp");

    // Thinning keeps the first, the latest and the pinned checkpoints, and seeking stays exact.
    const size_t before = timeline.checkpointCount();
    timeline.thin(timeline.bytes() / 2, ~size_t(0));
    bool kept = timeline.checkpointCount() < before && timeline.checkpointFor(0).step == 0 &&
                timeline.checkpointFor(kSteps).step > kSteps - 60 && timeline.checkpointFor(kAction).step == kAction;
    truth(kept, "timeline: thinning keeps the first, the latest and the pinned checkpoints");
    exact = true;
    for (size_t target : {size_t(99), size_t(222), size_t(371)}) {
        BodySystem b;
        exact = exact && seek(target, b) && sceneText(settings, camera, b) == expected[target];
    }
    truth(exact, "timeline: after thinning, still bit for bit");

    // A new future from the past: what came after is dropped, and recording goes on.
    {
        BodySystem b;
        seek(200, b);
        timeline.truncate();
        truth(timeline.atEnd() && timeline.length() == 200 && timeline.lastCheckpointTime() <= timeline.end(),
              "timeline: a new future drops the old one");
        stepScene(b, 0.5);
        timeline.addStep({0.5, 1.0f / 60.0f, {}}, b.time);
        truth(timeline.length() == 201 && timeline.end() == b.time, "timeline: and records on from there");
    }
}

// --- Camera keyframes (src/keyframes.h) -------------------------------------------------

static void testKeyframes() {
    using namespace keyframes;
    const double spin = 0.9;
    // Five shots around the hole, near and far, turned every way, with a hold at the third.
    std::vector<Keyframe> keys;
    const double place[5][7] = {
        // r, theta, phi, yaw, pitch, roll, fov
        {40.0, 1.40, 0.0, 0.0, 0.0, 0.0, 60.0},  {25.0, 1.10, 1.2, 0.3, -0.2, 0.1, 50.0},
        {12.0, 1.55, 2.6, -0.4, 0.1, -0.3, 70.0}, {30.0, 0.60, 4.0, 1.2, 0.5, 0.0, 40.0},
        {60.0, 1.30, 5.5, 0.0, -0.1, 0.4, 65.0},
    };
    bool cameraRoundTrip = true;
    for (int i = 0; i < 5; ++i) {
        Camera c;
        c.r = place[i][0];
        c.theta = place[i][1];
        c.phi = place[i][2];
        c.yaw = static_cast<float>(place[i][3]);
        c.pitch = static_cast<float>(place[i][4]);
        c.roll = static_cast<float>(place[i][5]);
        c.fovDeg = static_cast<float>(place[i][6]);
        Keyframe k = fromCamera(c, spin);
        k.duration = 3.0 + i;
        k.hold = i == 2 ? 1.5 : 0.0;
        keys.push_back(k);
        Pose pose;
        pose.position = k.position;
        pose.orientation = k.orientation;
        pose.fovDeg = k.fovDeg;
        Camera back;
        applyPose(pose, spin, 1.0, back);
        cameraRoundTrip = cameraRoundTrip && std::abs(back.r - c.r) < 1e-9 && std::abs(back.theta - c.theta) < 1e-9 &&
                          std::abs(std::remainder(back.phi - c.phi, 2.0 * kPi)) < 1e-9 &&
                          std::abs(back.yaw - c.yaw) < 1e-5f && std::abs(back.pitch - c.pitch) < 1e-5f &&
                          std::abs(back.roll - c.roll) < 1e-5f && back.fovDeg == c.fovDeg;
    }
    truth(cameraRoundTrip, "keyframes: a camera's place and view, as a keyframe and back");

    Path path;
    path.build(keys);
    // The spline passes through every keyframe, and the camera is there when due.
    double worst = 0.0;
    for (size_t i = 0; i < keys.size(); ++i) {
        worst = std::max(worst, glm::length(path.at(path.arrival(i)).position - keys[i].position));
        if (i + 1 < keys.size()) {
            worst = std::max(worst, glm::length(path.curve(i, 0.0) - keys[i].position));
            worst = std::max(worst, glm::length(path.curve(i, 1.0) - keys[i + 1].position));
        }
    }
    near(worst, 0.0, 1e-9, "keyframes: the spline passes through every keyframe");
    bool turned = true, zoomed = true;
    for (size_t i = 0; i < keys.size(); ++i) {
        const Pose p = path.at(path.arrival(i));
        turned = turned && std::abs(std::abs(glm::dot(p.orientation, keys[i].orientation)) - 1.0) < 1e-9;
        zoomed = zoomed && std::abs(p.fovDeg - keys[i].fovDeg) < 1e-9;
    }
    truth(turned && zoomed, "keyframes: the camera faces each keyframe's way, at its field of view");
    // Reparametrized by arc length, the camera moves evenly along each segment.
    double spread = 0.0;
    for (size_t i = 0; i < path.segments(); ++i) {
        const int n = 200;
        const double step = path.segmentLength(i) / n;
        for (int k = 0; k < n; ++k) {
            const double d = glm::length(path.evenly(i, (k + 1.0) / n) - path.evenly(i, static_cast<double>(k) / n));
            spread = std::max(spread, std::abs(d / step - 1.0));
        }
    }
    near(spread, 0.0, 0.01, "keyframes: arc-length speed varies by under 1%");
    // At rest at the start, through the hold and at the end; the distance never runs back.
    const Pose start = path.at(0.0), end = path.at(path.duration());
    const Pose held = path.at(path.arrival(2) + 0.7);
    bool still = glm::length(start.velocity) < 1e-9 && glm::length(end.velocity) < 1e-9 &&
                 glm::length(held.velocity) < 1e-9 && glm::length(held.position - keys[2].position) < 1e-9;
    truth(still, "keyframes: at rest at the start, while holding and at the end");
    bool onward = true;
    glm::dvec3 last = path.at(0.0).position;
    double travelled = 0.0;
    for (double t = 0.01; t <= path.duration(); t += 0.01) {
        const Pose p = path.at(t);
        travelled += glm::length(p.position - last);
        last = p.position;
        onward = onward && glm::dot(p.velocity, path.at(t + 1e-4).position - p.position) >= -1e-12;
    }
    double length = 0.0;
    for (size_t i = 0; i < path.segments(); ++i) length += path.segmentLength(i);
    truth(onward && std::abs(travelled / length - 1.0) < 1e-3, "keyframes: the camera never turns back along the path");
    // Squad turns smoothly through a keyframe: the rate of turning just before and just
    // after it is the same (C1), where linking slerps would jump.
    {
        const double t = path.arrival(1), h = 1e-4;
        const glm::dquat a = path.at(t - h).orientation, b = path.at(t).orientation, c = path.at(t + h).orientation;
        auto turn = [](const glm::dquat& x, const glm::dquat& y) {
            const glm::dquat d = glm::conjugate(x) * (glm::dot(x, y) < 0.0 ? -y : y);
            return glm::dvec3(d.x, d.y, d.z) * 2.0;
        };
        const glm::dvec3 before = turn(a, b) / h, after = turn(b, c) / h;
        near(glm::length(after - before) / std::max(glm::length(before), 1e-9), 0.0, 0.02,
             "keyframes: the camera turns smoothly through a keyframe (squad)");
    }
    // The ZAMO's own speed: a path that moves with the frame dragging is at rest in it.
    {
        bool capped = false;
        const glm::dvec3 p(8.0, 0.0, 0.0);
        const auto rtp = fromCartesian(spin, p.x, p.y, p.z);
        const double omega = metricAt(spin, rtp[0], rtp[1]).omega;
        const glm::dvec3 v = zamoVelocity(p, omega * glm::dvec3(-p.y, p.x, 0.0), spin, 0.99, capped);
        const glm::dvec3 fast = zamoVelocity(glm::dvec3(500.0, 0.0, 0.0), glm::dvec3(0.0, 3.0, 0.0), spin, 0.99, capped);
        truth(glm::length(v) < 1e-12 && capped && std::abs(glm::length(fast) - 0.99) < 1e-12,
              "keyframes: dragged with the ZAMO is at rest; faster than light is capped");
    }
    // Saved and read back.
    std::vector<Keyframe> back;
    const bool read = fromJson(toJson(keys), back) && back.size() == keys.size();
    bool same = read;
    for (size_t i = 0; same && i < keys.size(); ++i) {
        same = back[i].position == keys[i].position && back[i].orientation == keys[i].orientation &&
               back[i].fovDeg == keys[i].fovDeg && back[i].duration == keys[i].duration && back[i].hold == keys[i].hold;
    }
    truth(same, "keyframes: saved and read back bit for bit");
}

// --- What the Event Horizon Telescope sees (src/eht.h) --------------------------------------

static void testEht() {
    using namespace eht;
    const Target& m87 = targets()[0];
    const Target& sgra = targets()[1];
    // theta_g = G M / (c^2 D).
    near(angularScale(m87.massSolar, m87.distancePc), 3.82, 0.01, "EHT: M of M87* subtends 3.82 uas");
    near(angularScale(sgra.massSolar, sgra.distancePc), 5.13, 0.01, "EHT: M of Sgr A* subtends 5.13 uas");
    // The shadow of a hole without spin: 2 sqrt(27) M across.
    const std::vector<Point> circle = criticalCurve(0.0, 1.0, 3600);
    near(2.0 * criticalRadius(circle, 0.3), 2.0 * std::sqrt(27.0), 1e-4, "EHT: shadow 2 sqrt(27) M across at a = 0");
    near(2.0 * std::sqrt(27.0) * angularScale(m87.massSolar, m87.distancePc), 39.7, 0.1,
         "EHT: M87*'s shadow is 39.7 uas (the ring: 42 +- 3)");
    near(2.0 * std::sqrt(27.0) * angularScale(sgra.massSolar, sgra.distancePc), 53.3, 0.1,
         "EHT: Sgr A*'s shadow is 53.3 uas (the ring: 51.8 +- 2.3)");

    // The traced shadow's edge against the critical curve: a = 0.9 seen 17 degrees from the
    // axis, and a = 0 (bisecting along rays from the center of the image).
    auto edge = [](const View& cam, double psi) {
        double lo = 1.0, hi = 12.0;  // falls in at lo, escapes at hi
        for (int k = 0; k < 30; ++k) {
            const double mid = 0.5 * (lo + hi);
            const bool in = traceRay(cam, Flow{}, mid * std::cos(psi), mid * std::sin(psi), false).captured;
            (in ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    };
    double worst = 0.0;
    for (const auto& [a, deg] : {std::pair{0.9, 17.0}, std::pair{0.9, 60.0}, std::pair{0.0, 40.0}}) {
        View cam;
        cam.a = a;
        cam.inclination = deg * kPi / 180.0;
        const std::vector<Point> curve = criticalCurve(a, cam.inclination, 400);
        for (int k = 0; k < 8; ++k) {
            const double psi = 2.0 * kPi * (k + 0.25) / 8.0;
            const double want = criticalRadius(curve, psi);
            worst = std::max(worst, std::abs(edge(cam, psi) / want - 1.0));
        }
    }
    near(worst, 0.0, 2e-3, "EHT: the traced shadow follows the critical curve (a = 0.9 at 17 and 60 deg, a = 0)");

    // Synchrotron pieces.
    near(besselK2(1.0) / 1.6248388986351774, 1.0, 1e-5, "EHT: K2(1)");
    near(besselK2(0.5) / 7.5501835512408695, 1.0, 1e-5, "EHT: K2(0.5)");
    near(besselK2(10.0) / 2.1509817006932767e-05, 1.0, 1e-5, "EHT: K2(10)");
    near(besselK2(1.0 / 50.0) / (2.0 * 50.0 * 50.0), 1.0, 1e-3, "EHT: K2(1/theta) -> 2 theta^2 for hot electrons");
    near(brightnessTemperature(planck(230e9, 1e10), 230e9) / 1e10, 1.0, 1e-3,
         "EHT: a blackbody's brightness temperature is its temperature (Rayleigh-Jeans)");
    {
        // Doubling the density doubles the emissivity; the spectrum peaks near nu_s and
        // falls as exp(-X^1/3) far above it.
        const double j1 = synchrotronEmissivity(230e9, 1e5, 10.0, 10.0, 0.8);
        const double j2 = synchrotronEmissivity(230e9, 2e5, 10.0, 10.0, 0.8);
        const double far = synchrotronEmissivity(230e13, 1e5, 10.0, 10.0, 0.8);
        truth(std::abs(j2 / j1 - 2.0) < 1e-12 && far < 1e-10 * j1, "EHT: synchrotron emissivity scales with density, cuts off above nu_s");
    }

    // An image: calibrated to the observed flux, the approaching side brighter.
    {
        const int threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
        View cam;
        cam.a = 0.9;
        cam.inclination = 60.0 * kPi / 180.0;
        cam.size = 40;
        cam.massSolar = m87.massSolar;
        Flow flow;
        const double thetaG = angularScale(m87.massSolar, m87.distancePc);
        flow.n0 = calibrate(cam, flow, m87.flux, thetaG, threads);
        cam.size = 48;
        const Image image = render(cam, flow, threads);
        near(flux(image, thetaG) / m87.flux, 1.0, 0.08, "EHT: calibrated to M87*'s 0.5 Jy at 230 GHz");
        double left = 0.0, right = 0.0;
        for (int y = 0; y < image.size; ++y) {
            for (int x = 0; x < image.size; ++x) (x < image.size / 2 ? left : right) += image.intensity[static_cast<size_t>(y) * image.size + x];
        }
        // Bardeen's alpha (L = -alpha sin i): gas coming toward us is on the left.
        truth(left > 1.5 * right, "EHT: the side of the flow coming toward us is the brighter (Doppler beaming)");
        int shadowPixels = 0;
        for (unsigned char c : image.shadow) shadowPixels += c;
        const double shadowArea = shadowPixels * image.pixel() * image.pixel();
        const double curveArea = [&] {
            const std::vector<Point> curve = criticalCurve(cam.a, cam.inclination, 400);
            double area = 0.0;
            for (size_t k = 0; k < curve.size(); ++k) {
                const Point& p = curve[k];
                const Point& q = curve[(k + 1) % curve.size()];
                area += 0.5 * (p.alpha * q.beta - q.alpha * p.beta);
            }
            return std::abs(area);
        }();
        near(shadowArea / curveArea, 1.0, 0.05, "EHT: the image's shadow has the critical curve's area");
    }
}

// --- The Penrose process (src/penrose.h) ---------------------------------------------------------

static void testPenrose() {
    using namespace penrose;
    const double a = 0.99;
    const double rh = horizonRadius(a), re = ergosphereRadius(a, kPi / 2);
    // E = gamma (alpha + omega varpi v_phi) and L = gamma varpi v_phi agree with the ZAMO launch.
    {
        const double v[3] = {-0.3, 0.1, -0.5};
        Constants c;
        State st;
        launchFromZamo(a, 1.6, 1.4, 0.0, v, 1.0, c, st);
        const Zamo z = penrose::zamoAt(a, 1.6, 1.4);
        const Momentum p = momentum(1.0, v);
        near(energy(z, p), c.E, 1e-12, "Penrose: ZAMO-frame energy matches the geodesic's E");
        near(angularMomentum(z, p), c.L, 1e-12, "Penrose: ZAMO-frame angular momentum matches L");
    }
    // A body falling from rest far away (E0 = 1, L0 = 0), at r inside and just outside the ergosphere.
    auto fallingVelocity = [&](double r, double v0[3]) {
        Constants c{1.0, 0.0, 0.0, 1.0};
        State st{0.0, r, kPi / 2, 0.0, -std::sqrt(radialPotential(a, c, r)), 0.0};
        zamoVelocity(a, c, st, v0);
    };
    double v0[3];
    fallingVelocity(re * 1.05, v0);
    truth(!bestSplit(a, re * 1.05, kPi / 2, 1.0, v0, 0.95).valid,
          "Penrose: outside the ergosphere nothing can have negative energy");
    // At the closest approach of a parabolic orbit deep in the ergosphere, it can.
    const double rs = splitRadius(a);
    double Lp = 0.0;
    truth(parabolicPeriapsis(a, rs, Lp) && std::abs(radialPotential(a, Constants{1.0, Lp, 0.0, 1.0}, rs)) < 1e-9,
          "Penrose: the parabolic orbit turns round at the split radius (R = 0 there)");
    {
        const Constants c{1.0, Lp, 0.0, 1.0};
        const State st{0.0, rs, kPi / 2, 0.0, 0.0, 0.0};
        zamoVelocity(a, c, st, v0);
    }
    const Split best = bestSplit(a, rs, kPi / 2, 1.0, v0, 0.95);
    truth(best.valid && best.E1 < 0.0 && best.E2 > best.E0, "Penrose: inside, a split gives a negative-energy piece");
    near(best.m1 * best.E1 + best.m2 * best.E2, best.m0 * best.E0, 1e-12, "Penrose: energy is conserved at the split");
    near(best.m1 * best.L1 + best.m2 * best.L2, best.m0 * best.L0, 1e-12, "Penrose: angular momentum is conserved");
    near(best.efficiency, -best.m1 * best.E1 / (best.m0 * best.E0), 1e-12, "Penrose: the gain is -m1 E1");
    truth(best.efficiency > 0.0 && best.efficiency <= kMaxEfficiency, "Penrose: gains, but less than (sqrt 2 - 1) / 2");
    // Whatever falls in has E - Omega_H L >= 0 (the area theorem): random throws that plunge.
    {
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        const double omegaH = horizonOmega(1.0, a);
        double worst = 1.0;
        int plunged = 0;
        for (int k = 0; k < 4000; ++k) {
            double v[3] = {uni(rng), 0.3 * uni(rng), uni(rng)};
            const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            const double speed = 0.999 * std::abs(uni(rng));
            for (double& x : v) x *= speed / std::max(n, 1e-9);
            const double r = rh + (re - rh) * (0.5 + 0.5 * uni(rng));
            Constants c;
            State st;
            launchFromZamo(a, r, kPi / 2 + 0.3 * uni(rng), 0.0, v, 1.0, c, st);
            if (trajectory::classify(a, c, st).fate != trajectory::Fate::Plunges) continue;
            ++plunged;
            worst = std::min(worst, c.E - omegaH * c.L);
        }
        truth(plunged > 500 && worst >= -1e-12, "Penrose: everything that falls in has E - Omega_H L >= 0 (area theorem)");
    }
    // The hole: M^2 = M_irr^2 + J^2 / (4 M_irr^2), and 29% of a maximal hole is rotation.
    double worstIdentity = 0.0;
    for (double spin : {0.0, 0.3, 0.7, 0.9, 0.998}) {
        const double mirr = irreducibleMass(1.0, spin);
        worstIdentity = std::max(worstIdentity, std::abs(mirr * mirr + spin * spin / (4.0 * mirr * mirr) - 1.0));
    }
    near(worstIdentity, 0.0, 1e-12, "Penrose: Christodoulou's M^2 = M_irr^2 + J^2 / 4 M_irr^2");
    near(1.0 - irreducibleMass(1.0, 1.0), 1.0 - 1.0 / std::sqrt(2.0), 1e-12, "Penrose: 29.3% extractable at a = 1");
    near(horizonArea(1.0, 0.0), 16.0 * kPi, 1e-12, "Penrose: horizon area 16 pi M^2 without spin");
    // After the negative-energy piece falls in, the hole turns more slowly and its area grows.
    {
        const double after = spinAfter(a, best), mass = massAfter(best);
        truth(after < a && mass < 1.0 && horizonArea(mass, after * mass) >= horizonArea(1.0, a),
              "Penrose: the hole spins down, loses mass, and its area grows");
    }
    // The whole thing among the bodies: the body falls, splits, one piece falls in, the other leaves.
    {
        BodySystem system;
        system.addPreset(a, Preset::Penrose);
        bool escaped = false;
        for (int k = 0; k < 600 && !escaped; ++k) {
            system.update(a, 0.5, false);
            for (const Body& b : system.bodies) escaped = escaped || (b.name == "Escaping fragment" && b.state.r > 30.0);
        }
        const bool split = !system.penroseSplits.empty();
        truth(split && system.penroseSplits[0].absorbed && escaped && system.penroseSplits[0].E2 > system.penroseSplits[0].E0,
              "Penrose: the falling body splits; one piece falls in, the other escapes with more energy");
    }
}

// --- Naked singularities (spin above 1) -------------------------------------------------------

static void testNakedSingularity() {
    const double a = 1.2, b = std::sqrt(a * a - 1.0);
    // No horizon: Delta = (r - 1)^2 + b^2 > 0 everywhere, least b^2 at r = 1.
    double least = 1e30;
    for (double r = -20.0; r <= 20.0; r += 0.001) least = std::min(least, r * r - 2.0 * r + a * a);
    near(least, b * b, 1e-9, "naked: Delta > 0 everywhere, b^2 at its least");
    truth(horizonRadius(a) == 0.0 && innerHorizonRadius(a) == 0.0, "naked: no horizons");
    // The Kerr-Schild shifts in closed form have the right derivatives, for r < 0 too.
    double worst = 0.0;
    for (double r : {-6.0, -1.0, -0.2, 0.0, 0.3, 1.0, 2.5, 9.0}) {
        const double h = 1e-5, delta = r * r - 2.0 * r + a * a;
        const double dt = (timeShift(a, r + h) - timeShift(a, r - h)) / (2.0 * h);
        const double dphi = (azimuthShift(a, r + h) - azimuthShift(a, r - h)) / (2.0 * h);
        worst = std::max({worst, std::abs(dt - 2.0 * r / delta), std::abs(dphi - a / delta)});
    }
    near(worst, 0.0, 1e-8, "naked: d/dr of the closed-form shifts is 2r/Delta and a/Delta");
    // Equatorial photon orbits: r (r - 3)^2 = 4 a^2; past a = 1 only the retrograde one.
    const double rph = photonOrbitRadius(a, false);
    near(rph * (rph - 3.0) * (rph - 3.0), 4.0 * a * a, 1e-10, "naked: the retrograde photon orbit solves r (r-3)^2 = 4a^2");
    truth(photonOrbitRadius(a, true) == 0.0, "naked: no prograde equatorial photon orbit");
    // The ISCO: the least energy of circular orbits (dE/dr = 0).
    {
        const double r = iscoRadius(a, true), h = 1e-4;
        auto energy = [&](double x) { return circularOrbit(a, x, true).E; };
        near((energy(r + h) - energy(r - h)) / (2.0 * h), 0.0, 1e-6, "naked: dE/dr = 0 at the ISCO");
        truth(energy(r) < energy(r * 0.9) && energy(r) < energy(r * 1.1) && r < 1.0,
              "naked: the ISCO is a minimum of E, inside r = M");
    }
    // A ray through the ring's disk (r = 0 away from the equator) goes on to r < 0 on the
    // same geodesic: E, L, Q are kept (vr^2 = R, vtheta^2 = Theta all the way).
    {
        // Q < 0: such rays reach r = 0 (R(0) = -a^2 Q); with a^2 + Q small, R stays positive
        // for every r < 0 as well (otherwise the repulsion there turns the ray back).
        Constants c{1.0, 0.0, -1.2, 0.0};
        State st{0.0, 8.0, 0.3, 0.0, 0.0, 0.0};
        st.vr = -std::sqrt(radialPotential(a, c, st.r));
        st.vtheta = std::sqrt(std::max(polarPotential(a, c, st.theta), 0.0));
        double leastSigma = 1e30, drift = 0.0;
        for (int k = 0; k < 200000 && st.r > -8.0 && st.r < 1e3; ++k) {
            stepMino(a, c, st, photonStepSize(a, c, st, 0.005));
            enforceConstraints(a, c, st);
            normalizeAngles(st);
            const double cs = std::cos(st.theta);
            leastSigma = std::min(leastSigma, st.r * st.r + a * a * cs * cs);
            const double R = radialPotential(a, c, st.r), Th = polarPotential(a, c, st.theta);
            drift = std::max({drift, std::abs(st.vr * st.vr - R) / (1.0 + std::abs(R)),
                              std::abs(st.vtheta * st.vtheta - Th) / (1.0 + std::abs(Th))});
        }
        truth(st.r <= -8.0 && leastSigma > 0.05, "naked: a ray through the ring's disk comes out at r < 0");
        near(drift, 0.0, 1e-6, "naked: ... keeping E, L and Q (vr^2 = R, vtheta^2 = Theta)");
    }
}

static void run(const char* name, void (*test)()) {
    const auto start = std::chrono::steady_clock::now();
    test();
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (seconds > 1.0) std::printf("      (%s took %.1f s)\n", name, seconds);
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    run("testSceneFile", testSceneFile);
    run("testTimeline", testTimeline);
    run("testKeyframes", testKeyframes);
    run("testEht", testEht);
    run("testPenrose", testPenrose);
    run("testNakedSingularity", testNakedSingularity);
    run("testMerger", testMerger);
    run("testTidalFlare", testTidalFlare);
    run("testCentralUnits", testCentralUnits);
    run("testSurfaceAppearance", testSurfaceAppearance);
    run("testSurfacePhysics", testSurfacePhysics);
    run("testSurfaceLightBending", testSurfaceLightBending);
    run("testSurfaceRays", testSurfaceRays);
    run("testLens", testLens);
    run("testMercury", testMercury);
    run("testStellarStructure", testStellarStructure);
    run("testAccretionEfficiency", testAccretionEfficiency);
    run("testTidalScales", testTidalScales);
    run("testCollapse", testCollapse);
    run("testRealWorld", testRealWorld);
    run("testPulsarBeams", testPulsarBeams);
    run("testCloseBinary", testCloseBinary);
    run("testRocheLobe", testRocheLobe);
    run("testTidalModel", testTidalModel);
    run("testMutualGravity", testMutualGravity);
    run("testKerrSchildMetric", testKerrSchildMetric);
    run("testFallThroughHorizon", testFallThroughHorizon);
    run("testObserverFrames", testObserverFrames);
    run("testTrajectoryClassification", testTrajectoryClassification);
    run("testRayTrace", testRayTrace);
    run("testCartesianRoundTrip", testCartesianRoundTrip);
    run("testMovingObserver", testMovingObserver);
    run("testDebris", testDebris);
    run("testNearAxisPhotons", testNearAxisPhotons);
    run("testRadii", testRadii);
    run("testCircularOrbits", testCircularOrbits);
    run("testZamoRoundTrip", testZamoRoundTrip);
    run("testCircularIntegration", testCircularIntegration);
    run("testConstraintConservation", testConstraintConservation);
    run("testEquatorialOrbitSolver", testEquatorialOrbitSolver);
    run("testPeriapsisPrecession", testPeriapsisPrecession);
    run("testShadow", testShadow);
    std::printf("%d/%d checks passed\n", g_checks - g_failures, g_checks);
    return g_failures == 0 ? 0 : 1;
}
