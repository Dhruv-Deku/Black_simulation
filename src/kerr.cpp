#include "kerr.h"

#include <algorithm>
#include <cmath>

namespace kerr {

namespace {

// sin^3(theta) that keeps its sign but never reaches zero on the axis.
double safeSinCubed(double s) {
    double m = std::max(std::abs(s), 1e-9);
    return (s < 0.0 ? -1.0 : 1.0) * m * m * m;
}

struct Vec5 {
    double r, theta, phi, vr, vtheta;
};

Vec5 rateOf(const Derivatives& d, double scale) {
    return {d.dr * scale, d.dtheta * scale, d.dphi * scale, d.dvr * scale, d.dvtheta * scale};
}

State offset(const State& s, const Vec5& k, double h) {
    State o = s;
    o.r += h * k.r;
    o.theta += h * k.theta;
    o.phi += h * k.phi;
    o.vr += h * k.vr;
    o.vtheta += h * k.vtheta;
    return o;
}

// Equatorial (Q = 0, theta = pi/2) radial potential for a massive particle.
double equatorialR(double a, double E, double L, double r) {
    double P = E * (r * r + a * a) - a * L;
    double delta = r * r - 2.0 * r + a * a;
    double x = L - a * E;
    return P * P - delta * (r * r + x * x);
}

double equatorialHalfRPrime(double a, double E, double L, double r) {
    double P = E * (r * r + a * a) - a * L;
    double delta = r * r - 2.0 * r + a * a;
    double x = L - a * E;
    return 2.0 * E * r * P - (r - 1.0) * (r * r + x * x) - delta * r;
}

bool solveEquatorialOrbit(double a, double rp, double ra, bool prograde, double& E, double& L,
                          bool haveGuess) {
    if (!(rp > horizonRadius(a)) || !(ra > rp)) return false;
    const double sign = prograde ? 1.0 : -1.0;

    if (!haveGuess) {
        // Exact Schwarzschild solution as the starting point.
        double p = 2.0 * ra * rp / (ra + rp);
        double e = (ra - rp) / (ra + rp);
        double den = p - 3.0 - e * e;
        if (den > 1e-3) {
            E = std::sqrt(std::max(((p - 2.0) * (p - 2.0) - 4.0 * e * e) / (p * den), 1e-6));
            L = sign * p / std::sqrt(den);
        } else {
            E = 0.97;
            L = sign * 1.2 * std::sqrt(p);
        }
    }

    // Newton iteration on R(rp) = R(ra) = 0, each equation scaled by r^4.
    bool converged = false;
    for (int it = 0; it < 100; ++it) {
        double f[2], j[2][2];
        const double radii[2] = {rp, ra};
        for (int k = 0; k < 2; ++k) {
            double r = radii[k];
            double inv = 1.0 / (r * r * r * r);
            double P = E * (r * r + a * a) - a * L;
            double delta = r * r - 2.0 * r + a * a;
            double x = L - a * E;
            f[k] = equatorialR(a, E, L, r) * inv;
            j[k][0] = (2.0 * P * (r * r + a * a) + 2.0 * a * delta * x) * inv;  // dR/dE
            j[k][1] = (-2.0 * a * P - 2.0 * delta * x) * inv;                    // dR/dL
        }
        double det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if (std::abs(det) < 1e-300) return false;
        double dE = (f[0] * j[1][1] - f[1] * j[0][1]) / det;
        double dL = (j[0][0] * f[1] - j[1][0] * f[0]) / det;
        // Damp large steps so Newton cannot jump to a different branch.
        double damp = 1.0;
        damp = std::min(damp, 0.1 * std::max(std::abs(E), 0.1) / std::max(std::abs(dE), 1e-300));
        damp = std::min(damp, 0.1 * std::max(std::abs(L), 0.5) / std::max(std::abs(dL), 1e-300));
        E -= damp * dE;
        L -= damp * dL;
        if (std::abs(f[0]) < 1e-14 && std::abs(f[1]) < 1e-14 && damp == 1.0) {
            converged = true;
            break;
        }
    }
    if (!converged) return false;

    // A genuine bound orbit oscillating between rp and ra.
    if (!(E > 0.0 && E < 1.0) || L * sign <= 0.0) return false;
    if (equatorialHalfRPrime(a, E, L, rp) <= 0.0) return false;
    if (equatorialHalfRPrime(a, E, L, ra) >= 0.0) return false;
    for (int k = 1; k < 32; ++k) {
        double r = rp + (ra - rp) * k / 32.0;
        if (equatorialR(a, E, L, r) <= 0.0) return false;
    }
    return true;
}

}  // namespace

// --- Characteristic radii --------------------------------------------------

double horizonRadius(double a) {
    return 1.0 + std::sqrt(std::max(0.0, 1.0 - a * a));
}

double innerHorizonRadius(double a) {
    return 1.0 - std::sqrt(std::max(0.0, 1.0 - a * a));
}

double interiorLimit(double a) {
    return innerHorizonRadius(a) + 0.05;
}

double ergosphereRadius(double a, double theta) {
    double c = std::cos(theta);
    return 1.0 + std::sqrt(std::max(0.0, 1.0 - a * a * c * c));
}

double iscoRadius(double a, bool prograde) {
    double z1 = 1.0 + std::cbrt(1.0 - a * a) * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
    double z2 = std::sqrt(3.0 * a * a + z1 * z1);
    double root = std::sqrt(std::max(0.0, (3.0 - z1) * (3.0 + z1 + 2.0 * z2)));
    return prograde ? 3.0 + z2 - root : 3.0 + z2 + root;
}

double photonOrbitRadius(double a, bool prograde) {
    double s = prograde ? -1.0 : 1.0;
    return 2.0 * (1.0 + std::cos(2.0 / 3.0 * std::acos(s * a)));
}

double marginallyBoundRadius(double a, bool prograde) {
    double s = prograde ? 1.0 : -1.0;
    return 2.0 - s * a + 2.0 * std::sqrt(1.0 - s * a);
}

// --- Metric ------------------------------------------------------------------

Metric metricAt(double a, double r, double theta) {
    double s = std::sin(theta), c = std::cos(theta);
    Metric m;
    m.sigma = r * r + a * a * c * c;
    m.delta = r * r - 2.0 * r + a * a;
    double rra = r * r + a * a;
    m.A = rra * rra - a * a * m.delta * s * s;
    m.omega = 2.0 * a * r / m.A;
    return m;
}

MetricTensor metricTensor(double a, double r, double theta) {
    const double s = std::sin(theta), c = std::cos(theta);
    const double s2 = s * s;
    const double sigma = r * r + a * a * c * c;
    const double rra = r * r + a * a;
    const double A = rra * rra - a * a * (r * r - 2.0 * r + a * a) * s2;
    const double f = 2.0 * r / sigma;
    MetricTensor m{};
    m.g[0][0] = -(1.0 - f);
    m.g[0][1] = m.g[1][0] = f;
    m.g[0][3] = m.g[3][0] = -f * a * s2;
    m.g[1][1] = 1.0 + f;
    m.g[1][3] = m.g[3][1] = -a * s2 * (1.0 + f);
    m.g[2][2] = sigma;
    m.g[3][3] = A * s2 / sigma;
    return m;
}

double dot(const MetricTensor& g, const Vec4& x, const Vec4& y) {
    double sum = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) sum += g.g[i][j] * x[i] * y[j];
    }
    return sum;
}

Vec4 lower(const MetricTensor& g, const Vec4& x) {
    Vec4 out{};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) out[i] += g.g[i][j] * x[j];
    }
    return out;
}

double timeShift(double a, double r) {
    // integral of (r^2 + a^2) / Delta - 1, using r+^2 + a^2 = 2 r+ (and the same for r-).
    const double rp = horizonRadius(a), rm = innerHorizonRadius(a);
    const double inner = rm > 0.0 ? 2.0 * rm * std::log(std::abs(r - rm)) : 0.0;
    return (2.0 * rp * std::log(std::abs(r - rp)) - inner) / (rp - rm);
}

double azimuthShift(double a, double r) {
    const double rp = horizonRadius(a), rm = innerHorizonRadius(a);
    return a / (rp - rm) * std::log(std::abs((r - rp) / (r - rm)));
}

// --- Circular orbits -------------------------------------------------------------

double circularOmega(double a, double r, bool prograde) {
    double s = prograde ? 1.0 : -1.0;
    return s / (r * std::sqrt(r) + s * a);
}

double circularZamoSpeed(double a, double r, bool prograde) {
    double s = prograde ? 1.0 : -1.0;
    double sr = std::sqrt(r);
    double delta = r * r - 2.0 * r + a * a;
    return s * (r * r - 2.0 * s * a * sr + a * a) / (std::sqrt(delta) * (r * sr + s * a));
}

Constants circularOrbit(double a, double r, bool prograde) {
    double s = prograde ? 1.0 : -1.0;
    double sr = std::sqrt(r);
    double r32 = r * sr;
    double den = std::pow(r, 0.75) * std::sqrt(r32 - 3.0 * sr + 2.0 * s * a);
    Constants c;
    c.E = (r32 - 2.0 * sr + s * a) / den;
    c.L = s * (r * r - 2.0 * s * a * sr + a * a) / den;
    c.Q = 0.0;
    c.mu = 1.0;
    return c;
}

// --- Geodesic equations ------------------------------------------------------------

double inflowFactor(double P, double vr, double delta, double K) {
    if (P * vr < 0.0) return K / (P - vr);
    // Same signs: no cancellation. (Outward motion right at the horizon is genuinely
    // singular in these coordinates: it is the past horizon, which nothing crosses.)
    const double d = std::abs(delta) > 1e-12 ? delta : std::copysign(1e-12, delta);
    return (P + vr) / d;
}

Derivatives derivatives(double a, const Constants& c, const State& s) {
    const double r = s.r;
    const double sn = std::sin(s.theta), cs = std::cos(s.theta);
    const double sn2 = std::max(sn * sn, 1e-18);
    const double a2 = a * a, mu2 = c.mu * c.mu;
    const double delta = r * r - 2.0 * r + a2;
    const double rra = r * r + a2;
    const double P = c.E * rra - a * c.L;
    const double x = c.L - a * c.E;
    const double K = mu2 * r * r + x * x + c.Q;
    const double inflow = inflowFactor(P, s.vr, delta, K);

    Derivatives d;
    d.dr = s.vr;
    d.dtheta = s.vtheta;
    d.dvr = 2.0 * c.E * r * P - (r - 1.0) * K - delta * mu2 * r;
    d.dvtheta = sn * cs * a2 * (mu2 - c.E * c.E) + c.L * c.L * cs / safeSinCubed(sn);
    d.dphi = a * inflow - a * c.E + c.L / sn2;
    d.dt = rra * inflow - a * (a * c.E * sn * sn - c.L) - s.vr;
    return d;
}

double radialPotential(double a, const Constants& c, double r) {
    const double delta = r * r - 2.0 * r + a * a;
    const double P = c.E * (r * r + a * a) - a * c.L;
    const double x = c.L - a * c.E;
    return P * P - delta * (c.mu * c.mu * r * r + x * x + c.Q);
}

double polarPotential(double a, const Constants& c, double theta) {
    const double sn = std::sin(theta), cs = std::cos(theta);
    return c.Q - cs * cs * (a * a * (c.mu * c.mu - c.E * c.E) + c.L * c.L / std::max(sn * sn, 1e-18));
}

void launchFromZamo(double a, double r, double theta, double phi, const double v[3],
                    double mu, Constants& c, State& s) {
    const Metric m = metricAt(a, r, theta);
    const double sn = std::sin(theta), cs = std::cos(theta);

    // Local energy: gamma for a particle, 1 for a photon of unit ZAMO frequency.
    double gamma = 1.0;
    if (mu > 0.0) {
        double v2 = std::min(v[0] * v[0] + v[1] * v[1] + v[2] * v[2], 1.0 - 1e-12);
        gamma = 1.0 / std::sqrt(1.0 - v2);
    }

    c.mu = mu;
    c.L = gamma * v[2] * std::sqrt(m.A / m.sigma) * sn;
    c.E = gamma * std::sqrt(m.sigma * m.delta / m.A) + m.omega * c.L;

    s.t = 0.0;
    s.r = r;
    s.theta = theta;
    s.phi = phi;
    s.vr = gamma * v[0] * std::sqrt(m.sigma * m.delta);  // = delta * p_r
    s.vtheta = gamma * v[1] * std::sqrt(m.sigma);         // = p_theta

    if (mu == 0.0) {
        // The photon's overall scale is arbitrary; fix E = 1.
        double inv = 1.0 / c.E;
        c.E = 1.0;
        c.L *= inv;
        s.vr *= inv;
        s.vtheta *= inv;
    }
    c.Q = s.vtheta * s.vtheta +
          cs * cs * (a * a * (mu * mu - c.E * c.E) + c.L * c.L / std::max(sn * sn, 1e-18));

    if (mu == 0.0 && std::abs(c.L) < kAxisSnap * std::sqrt(std::max(c.Q + a * a, 0.0))) {
        // This photon would pass within ~kAxisSnap radians of the spin axis, where
        // phi swings through pi faster than any step can resolve. Through the axis
        // itself (L = 0) the motion is smooth, so send it through the axis instead;
        // the path moves by a negligible amount.
        c.L = 0.0;
        c.Q = s.vtheta * s.vtheta - cs * cs * a * a;
    }
}

void zamoVelocity(double a, const Constants& c, const State& s, double v[3]) {
    const Metric m = metricAt(a, s.r, s.theta);
    const double sn = std::sin(s.theta);
    const double gamma = std::sqrt(m.A / (m.sigma * m.delta)) * (c.E - m.omega * c.L);
    v[0] = s.vr / std::sqrt(m.sigma * m.delta) / gamma;
    v[1] = s.vtheta / std::sqrt(m.sigma) / gamma;
    v[2] = c.L / (std::sqrt(m.A / m.sigma) * sn * gamma);
}

Vec4 fourVelocity(double a, const Constants& c, const State& s) {
    const Derivatives d = derivatives(a, c, s);
    const double invSigma = 1.0 / metricAt(a, s.r, s.theta).sigma;
    return {d.dt * invSigma, d.dr * invSigma, d.dtheta * invSigma, d.dphi * invSigma};
}

Vec4 zamoFourVelocity(double a, double r, double theta) {
    // u = (e_t + omega e_phi) / alpha. At fixed r the Kerr-Schild t and phi tick
    // like the Boyer-Lindquist ones, so the components are the same.
    const Metric m = metricAt(a, r, theta);
    const double ut = std::sqrt(m.A / (m.sigma * m.delta));
    return {ut, 0.0, 0.0, m.omega * ut};
}

Frame observerFrame(double a, double r, double theta, const Vec4& u) {
    // Gram-Schmidt from the observer's 4-velocity, then the theta, phi and r
    // directions. For a ZAMO, whose angular momentum u_phi is zero, this gives
    // exactly the ZAMO frame; for a radial or azimuthal boost, the boosted frame.
    const MetricTensor g = metricTensor(a, r, theta);
    Frame f;
    f.e[0] = u;
    const Vec4 candidates[3] = {{0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 0.0}};
    const int slot[3] = {2, 3, 1};
    for (int k = 0; k < 3; ++k) {
        Vec4 w = candidates[k];
        // Remove the parts along u (timelike, norm -1) and the axes found so far.
        const double along = dot(g, w, u);
        for (int i = 0; i < 4; ++i) w[i] += along * u[i];
        for (int j = 0; j < k; ++j) {
            const Vec4& e = f.e[slot[j]];
            const double p = dot(g, w, e);
            for (int i = 0; i < 4; ++i) w[i] -= p * e[i];
        }
        const double n = std::sqrt(std::max(dot(g, w, w), 1e-300));
        for (double& x : w) x /= n;
        f.e[slot[k]] = w;
    }
    return f;
}

void photonFromFrame(double a, double r, double theta, double phi, const Frame& f, const double look[3],
                     Constants& c, State& s) {
    // The photon arrives travelling along -look, with unit energy in this frame.
    Vec4 p;
    for (int i = 0; i < 4; ++i) {
        p[i] = f.e[0][i] - look[0] * f.e[1][i] - look[1] * f.e[2][i] - look[2] * f.e[3][i];
    }
    const Vec4 low = lower(metricTensor(a, r, theta), p);
    const double sigma = metricAt(a, r, theta).sigma;
    const double sn = std::sin(theta), cs = std::cos(theta);
    c.mu = 0.0;
    c.E = -low[0];
    c.L = low[3];
    s.t = 0.0;
    s.r = r;
    s.theta = theta;
    s.phi = phi;
    s.vr = sigma * p[1];
    s.vtheta = low[2];  // = Sigma p^theta
    c.Q = s.vtheta * s.vtheta + cs * cs * (c.L * c.L / std::max(sn * sn, 1e-18) - a * a * c.E * c.E);
    if (std::abs(c.L) < kAxisSnap * std::sqrt(std::max(c.Q + a * a * c.E * c.E, 0.0))) {
        // See launchFromZamo: through the axis instead of within a hair of it.
        c.L = 0.0;
        c.Q = s.vtheta * s.vtheta - cs * cs * a * a * c.E * c.E;
    }
}

double observedPhoton(double a, double r, double theta, double phi, const double vObs[3], const double look[3],
                      Constants& c, State& s) {
    Constants oc;
    State os;
    launchFromZamo(a, r, theta, phi, vObs, 1.0, oc, os);
    const Frame f = observerFrame(a, r, theta, fourVelocity(a, oc, os));
    photonFromFrame(a, r, theta, phi, f, look, c, s);
    // Normalize to E = 1; the observer then measures 1 / E_old.
    const double inv = 1.0 / c.E;
    c.L *= inv;
    c.Q *= inv * inv;
    s.vr *= inv;
    s.vtheta *= inv;
    c.E = 1.0;
    return inv;
}

void stepMino(double a, const Constants& c, State& s, double h) {
    auto f = [&](const State& x) { return derivatives(a, c, x); };
    const Derivatives d1 = f(s);
    const Vec5 k1 = rateOf(d1, 1.0);
    const Derivatives d2 = f(offset(s, k1, 0.5 * h));
    const Vec5 k2 = rateOf(d2, 1.0);
    const Derivatives d3 = f(offset(s, k2, 0.5 * h));
    const Vec5 k3 = rateOf(d3, 1.0);
    const Derivatives d4 = f(offset(s, k3, h));
    const Vec5 k4 = rateOf(d4, 1.0);

    s.t += h / 6.0 * (d1.dt + 2.0 * d2.dt + 2.0 * d3.dt + d4.dt);
    s.r += h / 6.0 * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
    s.theta += h / 6.0 * (k1.theta + 2.0 * k2.theta + 2.0 * k3.theta + k4.theta);
    s.phi += h / 6.0 * (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi);
    s.vr += h / 6.0 * (k1.vr + 2.0 * k2.vr + 2.0 * k3.vr + k4.vr);
    s.vtheta += h / 6.0 * (k1.vtheta + 2.0 * k2.vtheta + 2.0 * k3.vtheta + k4.vtheta);
}

void enforceConstraints(double a, const Constants& c, State& s) {
    // Fraction of the potential's scale below which we are "near a turning point".
    constexpr double kTurningZone = 0.05;

    const double r = s.r;
    const double delta = r * r - 2.0 * r + a * a;
    const double P = c.E * (r * r + a * a) - a * c.L;
    const double x = c.L - a * c.E;
    const double K = c.mu * c.mu * r * r + x * x + c.Q;
    const double R = P * P - delta * K;
    if (R > kTurningZone * (P * P + std::abs(delta) * K)) {
        s.vr = std::copysign(std::sqrt(R), s.vr);
    }

    const double sn = std::sin(s.theta), cs = std::cos(s.theta);
    const double cot2 = cs * cs / std::max(sn * sn, 1e-18);
    const double polarMass = a * a * (c.mu * c.mu - c.E * c.E) * cs * cs;
    const double Theta = c.Q - polarMass - c.L * c.L * cot2;
    if (Theta > kTurningZone * (std::abs(c.Q) + std::abs(polarMass) + c.L * c.L * cot2)) {
        s.vtheta = std::copysign(std::sqrt(Theta), s.vtheta);
    }
}

void stepCoordinateTime(double a, const Constants& c, State& s, double h) {
    // Same system with t as the independent variable: dy/dt = (dy/dtau) / (dt/dtau).
    auto f = [&](const State& x) {
        Derivatives d = derivatives(a, c, x);
        return rateOf(d, 1.0 / d.dt);
    };
    const Vec5 k1 = f(s);
    const Vec5 k2 = f(offset(s, k1, 0.5 * h));
    const Vec5 k3 = f(offset(s, k2, 0.5 * h));
    const Vec5 k4 = f(offset(s, k3, h));

    s.t += h;
    s.r += h / 6.0 * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
    s.theta += h / 6.0 * (k1.theta + 2.0 * k2.theta + 2.0 * k3.theta + k4.theta);
    s.phi += h / 6.0 * (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi);
    s.vr += h / 6.0 * (k1.vr + 2.0 * k2.vr + 2.0 * k3.vr + k4.vr);
    s.vtheta += h / 6.0 * (k1.vtheta + 2.0 * k2.vtheta + 2.0 * k3.vtheta + k4.vtheta);
}

double photonStepSize(double a, const Constants& c, const State& s, double tolerance) {
    const Derivatives d = derivatives(a, c, s);
    const double sn = std::abs(std::sin(s.theta));
    const double tol = tolerance * std::clamp(s.r / 8.0, 1.0, 4.0);
    double h = tol / (std::abs(d.dr) / s.r + std::abs(d.dtheta) / std::max(sn, 1e-4) + std::abs(d.dphi) + 1e-6);
    h = std::min(h, std::sqrt(tol * s.r / (std::abs(d.dvr) + 1e-6)));
    h = std::min(h, std::sqrt(tol * std::max(sn, 1e-4) / (std::abs(d.dvtheta) + 1e-6)));
    return h;
}

double suggestedTimeStep(double a, const Constants& c, const State& s, double tolerance) {
    const Derivatives d = derivatives(a, c, s);
    const double inv = 1.0 / d.dt;
    const double sn = std::max(std::abs(std::sin(s.theta)), 1e-3);

    // Velocity limit: relative change of r, sin(theta) and phi per step. Steps
    // shrink approaching the spin axis, where phi can swing through ~pi quickly.
    double rate = std::abs(d.dr * inv) / s.r + std::abs(d.dtheta * inv) / sn + std::abs(d.dphi * inv);
    double h = tolerance / std::max(rate, 1e-12);

    // Acceleration limit, for particles momentarily at rest.
    double acc = std::abs(d.dvr) * inv * inv / s.r + std::abs(d.dvtheta) * inv * inv / sn;
    h = std::min(h, std::sqrt(tolerance / std::max(acc, 1e-30)));

    // Never more than a fraction of the local orbital period.
    return std::min(h, 0.05 * 2.0 * kPi * s.r * std::sqrt(s.r));
}

void normalizeAngles(State& s) {
    s.theta = std::fmod(s.theta, 2.0 * kPi);
    if (s.theta < 0.0) s.theta += 2.0 * kPi;
    if (s.theta > kPi) {
        // Passing over the pole: (theta, phi) -> (2 pi - theta, phi + pi).
        s.theta = 2.0 * kPi - s.theta;
        s.phi += kPi;
        s.vtheta = -s.vtheta;
    }
}

// --- Equatorial bound orbits -----------------------------------------------------

bool equatorialOrbit(double a, double rp, double ra, bool prograde, double& E, double& L) {
    return solveEquatorialOrbit(a, rp, ra, prograde, E, L, false);
}

bool separatrixPeriapsis(double a, double ra, bool prograde, double& rp) {
    double hi = iscoRadius(a, prograde) * 1.001;
    double lo = photonOrbitRadius(a, prograde);
    if (ra <= hi * 1.01) return false;

    double E = 0.0, L = 0.0;
    if (!solveEquatorialOrbit(a, hi, ra, prograde, E, L, false)) return false;

    // Bisect on validity, continuing from the last valid (E, L).
    for (int it = 0; it < 60; ++it) {
        double mid = 0.5 * (lo + hi);
        double e = E, l = L;
        if (solveEquatorialOrbit(a, mid, ra, prograde, e, l, true)) {
            hi = mid;
            E = e;
            L = l;
        } else {
            lo = mid;
        }
    }
    rp = hi;
    return true;
}

// --- Cartesian coordinates ---------------------------------------------------------

std::array<double, 3> toCartesian(double a, double r, double theta, double phi) {
    const double sn = std::sin(theta), cp = std::cos(phi), sp = std::sin(phi);
    return {(r * cp - a * sp) * sn, (r * sp + a * cp) * sn, r * std::cos(theta)};
}

std::array<double, 3> fromCartesian(double a, double x, double y, double z) {
    // r solves (x^2 + y^2) / (r^2 + a^2) + z^2 / r^2 = 1.
    const double w = x * x + y * y + z * z - a * a;
    const double r = std::sqrt(0.5 * (w + std::sqrt(w * w + 4.0 * a * a * z * z)));
    const double theta = std::acos(std::clamp(r > 0.0 ? z / r : 1.0, -1.0, 1.0));
    return {r, theta, std::atan2(y, x) - std::atan2(a, r)};
}

std::array<double, 3> cartesianVelocity(double a, double r, double theta, double phi, double dr, double dtheta,
                                        double dphi) {
    // d(x + i y) = [dr sin(theta) + (r + i a)(cos(theta) dtheta + i sin(theta) dphi)] e^{i phi}
    const double sn = std::sin(theta), cs = std::cos(theta), cp = std::cos(phi), sp = std::sin(phi);
    const double re = dr * sn + r * cs * dtheta - a * sn * dphi;
    const double im = a * cs * dtheta + r * sn * dphi;
    return {re * cp - im * sp, re * sp + im * cp, dr * cs - r * sn * dtheta};
}

LocalAxes localAxes(double a, double r, double theta, double phi) {
    // Spherical axes at the point's Cartesian azimuth.
    const double az = phi + std::atan2(a, r);
    const double st = std::sin(theta), ct = std::cos(theta), sp = std::sin(az), cp = std::cos(az);
    return {{st * cp, st * sp, ct}, {ct * cp, ct * sp, -st}, {-sp, cp, 0.0}};
}

}  // namespace kerr
