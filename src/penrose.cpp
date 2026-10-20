#include "penrose.h"

#include "kerr.h"
#include "trajectory.h"

#include <cmath>

namespace penrose {

Zamo zamoAt(double a, double r, double theta) {
    const kerr::Metric m = kerr::metricAt(a, r, theta);
    return {std::sqrt(m.sigma * m.delta / m.A), m.omega, std::sqrt(m.A / m.sigma) * std::sin(theta)};
}

double Momentum::mass() const {
    const double m2 = p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3];
    return m2 >= 0.0 ? std::sqrt(m2) : -1.0;
}

Momentum momentum(double mass, const double v[3]) {
    const double gamma = 1.0 / std::sqrt(1.0 - (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
    Momentum p;
    p.p[0] = mass * gamma;
    for (int i = 0; i < 3; ++i) p.p[i + 1] = mass * gamma * v[i];
    return p;
}

double energy(const Zamo& z, const Momentum& p) {
    return z.alpha * p.p[0] + z.omega * z.varpi * p.p[3];
}

double angularMomentum(const Zamo& z, const Momentum& p) {
    return z.varpi * p.p[3];
}

Split split(double a, double r, double theta, double m0, const double v0[3], double m1, const double v1[3]) {
    Split s;
    s.r = r;
    s.theta = theta;
    s.spinBefore = a;
    s.m0 = m0;
    s.m1 = m1;
    const Zamo z = zamoAt(a, r, theta);
    const Momentum p0 = momentum(m0, v0), p1 = momentum(m1, v1);
    Momentum p2;
    for (int i = 0; i < 4; ++i) p2.p[i] = p0.p[i] - p1.p[i];
    s.m2 = p2.mass();
    if (s.m2 <= 0.0 || p2.p[0] <= 0.0) {
        s.why = "what is left would move faster than light";
        return s;
    }
    s.E0 = energy(z, p0) / m0;
    s.E1 = energy(z, p1) / m1;
    s.E2 = energy(z, p2) / s.m2;
    s.L0 = angularMomentum(z, p0) / m0;
    s.L1 = angularMomentum(z, p1) / m1;
    s.L2 = angularMomentum(z, p2) / s.m2;
    for (int i = 0; i < 3; ++i) {
        s.v1[i] = v1[i];
        s.v2[i] = p2.p[i + 1] / p2.p[0];
    }
    s.efficiency = (s.m2 * s.E2 - m0 * s.E0) / (m0 * s.E0);
    s.valid = true;
    return s;
}

Split bestSplit(double a, double r, double theta, double m0, const double v0[3], double speed) {
    Split best;
    best.why = "no split here gives a fragment of negative energy that falls in while the other escapes";
    const Zamo z = zamoAt(a, r, theta);
    if (z.alpha >= z.omega * z.varpi) {
        best.why = "outside the ergosphere: nothing can have negative energy here";
        return best;
    }
    auto fate = [&](const double v[3]) {
        kerr::Constants c;
        kerr::State st;
        kerr::launchFromZamo(a, r, theta, 0.0, v, 1.0, c, st);
        return trajectory::classify(a, c, st).fate;
    };
    // Fragment 1 thrown backwards (against the rotation), at angles from straight back
    // toward the hole; its share of the rest mass from a few percent to most.
    for (int k = 0; k <= 60; ++k) {
        const double angle = -kerr::kPi / 2.0 + kerr::kPi * k / 60.0;  // from inward to outward
        const double v1[3] = {speed * std::sin(angle), 0.0, -speed * std::cos(angle)};
        for (int j = 1; j < 40; ++j) {
            const double m1 = m0 * j / 40.0;
            Split s = split(a, r, theta, m0, v0, m1, v1);
            if (!s.valid || s.E1 >= 0.0 || s.efficiency <= best.efficiency) continue;
            if (fate(v1) != trajectory::Fate::Plunges || fate(s.v2) != trajectory::Fate::Escapes) continue;
            best = s;
        }
    }
    return best;
}

bool parabolicPeriapsis(double a, double r, double& L) {
    const double A = 2.0 - r, B = -4.0 * a, C = 2.0 * (r * r + a * a);
    const double D = B * B - 4.0 * A * C;
    if (D < 0.0 || std::abs(A) < 1e-12) return false;
    L = (-B - std::sqrt(D)) / (2.0 * A);
    return true;
}

double splitRadius(double a) {
    const double rh = kerr::horizonRadius(a), re = kerr::ergosphereRadius(a, 0.5 * kerr::kPi);
    return std::max(re - 0.7 * (re - rh), 1.02 * kerr::marginallyBoundRadius(a, true));
}

double horizonArea(double M, double a) {
    const double rp = M + std::sqrt(std::max(M * M - a * a, 0.0));
    return 8.0 * kerr::kPi * M * rp;
}

double irreducibleMass(double M, double a) {
    return std::sqrt(horizonArea(M, a) / (16.0 * kerr::kPi));
}

double horizonOmega(double M, double a) {
    const double rp = M + std::sqrt(std::max(M * M - a * a, 0.0));
    return a / (2.0 * M * rp);
}

double massAfter(const Split& s) {
    return 1.0 + s.m1 * s.E1;
}

double spinAfter(double a, const Split& s) {
    const double M = massAfter(s);
    return (a + s.m1 * s.L1) / (M * M);
}

}  // namespace penrose
