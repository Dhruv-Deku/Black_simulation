#include "debris.h"

#include <algorithm>
#include <cmath>

namespace debris {

namespace {

using Vec3 = std::array<double, 3>;

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Initial conditions are built in Boyer-Lindquist Cartesian coordinates,
// x + i y = sqrt(r^2 + a^2) e^{i phi_BL} sin(theta), whose spherical axes are exactly
// the ZAMO's azimuthal directions; positions are converted at the end.
Vec3 blCartesian(double a, double r, double theta, double phiBL) {
    const double rho = std::sqrt(r * r + a * a), sn = std::sin(theta);
    return {rho * sn * std::cos(phiBL), rho * sn * std::sin(phiBL), r * std::cos(theta)};
}

kerr::LocalAxes blAxes(double theta, double phiBL) {
    const double st = std::sin(theta), ct = std::cos(theta), sp = std::sin(phiBL), cp = std::cos(phiBL);
    return {{st * cp, st * sp, ct}, {ct * cp, ct * sp, -st}, {-sp, cp, 0.0}};
}

// Launches one particle at Boyer-Lindquist Cartesian position p with Cartesian
// velocity vel, given as components along the local axes. Returns false if it would
// start inside or at the horizon, or faster than light.
bool launch(double a, const Vec3& p, const Vec3& vel, float temperature, ParticleInit& out) {
    // Same r and theta as the Kerr-Schild mapping; only the azimuth differs.
    const double w = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - a * a;
    const double r = std::sqrt(0.5 * (w + std::sqrt(w * w + 4.0 * a * a * p[2] * p[2])));
    if (r < kerr::horizonRadius(a) + 0.1) return false;
    const double theta = std::acos(std::clamp(p[2] / r, -1.0, 1.0));
    const double phiBL = std::atan2(p[1], p[0]);
    const kerr::LocalAxes ax = blAxes(theta, phiBL);
    const double v[3] = {dot(vel, ax.r), dot(vel, ax.theta), dot(vel, ax.phi)};
    if (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] >= 0.999) return false;

    kerr::Constants c;
    kerr::State s;
    kerr::launchFromZamo(a, r, theta, phiBL + kerr::azimuthShift(a, r), v, 1.0, c, s);
    out = {static_cast<float>(s.r), static_cast<float>(s.theta), static_cast<float>(s.phi),
           static_cast<float>(s.vr), static_cast<float>(s.vtheta), static_cast<float>(c.E),
           static_cast<float>(c.L), static_cast<float>(c.Q), temperature};
    return true;
}

}  // namespace

std::vector<ParticleInit> cloud(double a, double r, double theta, double phi, const double v[3], double radius,
                                bool gaussian, int count, float temperature, std::mt19937& rng) {
    const double phiBL = phi - kerr::azimuthShift(a, r);
    const Vec3 center = blCartesian(a, r, theta, phiBL);
    const kerr::LocalAxes ax = blAxes(theta, phiBL);
    Vec3 vel;
    for (int i = 0; i < 3; ++i) vel[i] = v[0] * ax.r[i] + v[1] * ax.theta[i] + v[2] * ax.phi[i];

    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    std::uniform_real_distribution<float> warmth(0.85f, 1.15f);
    std::vector<ParticleInit> out;
    out.reserve(count);
    for (int attempt = 0; static_cast<int>(out.size()) < count && attempt < count * 4; ++attempt) {
        Vec3 offset;
        if (gaussian) {
            for (double& o : offset) o = normal(rng) * radius;
        } else {
            do {
                for (double& o : offset) o = uni(rng);
            } while (dot(offset, offset) > 1.0);
            for (double& o : offset) o *= radius;
        }
        const Vec3 p = {center[0] + offset[0], center[1] + offset[1], center[2] + offset[2]};
        ParticleInit init;
        if (launch(a, p, vel, temperature * warmth(rng), init)) out.push_back(init);
    }
    return out;
}

std::vector<ParticleInit> deformedStar(double a, double r, double theta, double phi, const double v[3],
                                       double radius, const tidal::Mat3& F, const tidal::Mat3& Fdot, int count,
                                       float temperature, std::mt19937& rng) {
    const double phiBL = phi - kerr::azimuthShift(a, r);
    const Vec3 center = blCartesian(a, r, theta, phiBL);
    const kerr::LocalAxes ax = blAxes(theta, phiBL);
    Vec3 vel;
    for (int i = 0; i < 3; ++i) vel[i] = v[0] * ax.r[i] + v[1] * ax.theta[i] + v[2] * ax.phi[i];

    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    std::uniform_real_distribution<float> warmth(0.9f, 1.1f);
    std::vector<ParticleInit> out;
    out.reserve(count);
    for (int attempt = 0; static_cast<int>(out.size()) < count && attempt < count * 4; ++attempt) {
        Vec3 X;
        do {
            for (double& o : X) o = uni(rng);
        } while (dot(X, X) > 1.0);
        for (double& o : X) o *= radius;
        Vec3 p, u;
        for (int i = 0; i < 3; ++i) {
            p[i] = center[i] + F[i][0] * X[0] + F[i][1] * X[1] + F[i][2] * X[2];
            u[i] = vel[i] + Fdot[i][0] * X[0] + Fdot[i][1] * X[1] + Fdot[i][2] * X[2];
        }
        ParticleInit init;
        if (launch(a, p, u, temperature * warmth(rng), init)) out.push_back(init);
    }
    return out;
}

std::vector<ParticleInit> ring(double a, double rInner, double rOuter, double tiltDeg, int count,
                               float temperature, std::mt19937& rng) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> jitter(0.0, 1.0);
    std::uniform_real_distribution<float> warmth(0.8f, 1.2f);
    const double tilt = tiltDeg * kerr::kPi / 180.0;
    const double ct = std::cos(tilt), st = std::sin(tilt);
    const bool prograde = std::abs(tiltDeg) <= 90.0;

    std::vector<ParticleInit> out;
    out.reserve(count);
    for (int attempt = 0; static_cast<int>(out.size()) < count && attempt < count * 4; ++attempt) {
        // Uniform in area between the two radii.
        const double R = std::sqrt(rInner * rInner + (rOuter * rOuter - rInner * rInner) * uni(rng));
        const double psi = 2.0 * kerr::kPi * uni(rng);
        const double c = std::cos(psi), s = std::sin(psi);
        // Circle of radius R in the plane tilted about the x axis, slightly thickened.
        const double height = 0.01 * R * jitter(rng);
        const Vec3 p = {R * c, R * s * ct - height * st, R * s * st + height * ct};
        const double speed = std::abs(kerr::circularZamoSpeed(a, R, prograde)) * (1.0 + 0.003 * jitter(rng));
        const Vec3 vel = {-s * speed, c * ct * speed, c * st * speed};
        ParticleInit init;
        if (launch(a, p, vel, temperature * warmth(rng), init)) out.push_back(init);
    }
    return out;
}

}  // namespace debris
