#include "gravity.h"

#include "central.h"

#include <algorithm>
#include <cmath>

namespace gravity {

namespace {

// Bodies lighter than this barely pull on anything (a star around a million-solar-mass
// hole is about 1e-6); they are still pulled by heavier ones.
constexpr double kNegligibleMass = 1e-5;

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>;

// d(x, y, z) / d(r, theta, phi) of the Kerr-Schild Cartesian mapping (columns: r, theta, phi).
Mat3 jacobian(double a, double r, double theta, double phi) {
    const double sn = std::sin(theta), cs = std::cos(theta), sp = std::sin(phi), cp = std::cos(phi);
    return {{{cp * sn, (r * cp - a * sp) * cs, -(r * sp + a * cp) * sn},
             {sp * sn, (r * sp + a * cp) * cs, (r * cp - a * sp) * sn},
             {cs, -r * sn, 0.0}}};
}

Vec3 solve(const Mat3& m, const Vec3& b) {
    const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                       m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                       m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    Vec3 x;
    for (int c = 0; c < 3; ++c) {
        Mat3 t = m;
        for (int r = 0; r < 3; ++r) t[r][c] = b[r];
        x[c] = (t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) - t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0]) +
                t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0])) / det;
    }
    return x;
}

}  // namespace

kerr::Vec4 pull(double a, const std::vector<Member>& members, size_t i) {
    const Member& mi = members[i];
    const kerr::State& si = mi.s;
    const Vec3 xi = kerr::toCartesian(a, si.r, si.theta, si.phi);
    const Mat3 J = jacobian(a, si.r, si.theta, si.phi);
    const kerr::Vec4 u = kerr::fourVelocity(a, mi.c, si);
    const kerr::MetricTensor g = kerr::metricTensor(a, si.r, si.theta);
    const kerr::Vec4 uLow = kerr::lower(g, u);
    kerr::Vec4 force{};
    for (size_t j = 0; j < members.size(); ++j) {
        const Member& mj = members[j];
        if (j == i || mj.mass <= 0.0) continue;
        const Vec3 xj = kerr::toCartesian(a, mj.s.r, mj.s.theta, mj.s.phi);
        const Vec3 dq = solve(J, {xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2]});
        kerr::Vec4 sep = {0.0, dq[0], dq[1], dq[2]};
        // Where body j is at the moment simultaneous with body i in i's own rest frame:
        // slide along j's worldline until the separation is orthogonal to u_i.
        const kerr::Derivatives dj = kerr::derivatives(a, mj.c, mj.s);
        const kerr::Vec4 vj = {1.0, dj.dr / dj.dt, dj.dtheta / dj.dt, dj.dphi / dj.dt};
        double us = 0.0, uv = 0.0;
        for (int k = 0; k < 4; ++k) {
            us += uLow[k] * sep[k];
            uv += uLow[k] * vj[k];
        }
        const double shift = -us / uv;
        for (int k = 0; k < 4; ++k) sep[k] += vj[k] * shift;
        const kerr::Vec4 sepLow = kerr::lower(g, sep);
        const double eps = 0.5 * (mi.softening + mj.softening);
        double d2 = eps * eps;
        for (int k = 0; k < 4; ++k) d2 += sepLow[k] * sep[k];
        const double f = mj.mass / (d2 * std::sqrt(d2));
        for (int k = 0; k < 4; ++k) force[k] += f * sepLow[k];
    }
    return force;
}

void kick(double a, std::vector<Member>& members, double dt) {
    // Measure everything before changing anything.
    std::vector<kerr::Vec4> force(members.size());
    for (size_t i = 0; i < members.size(); ++i) force[i] = pull(a, members, i);

    for (size_t i = 0; i < members.size(); ++i) {
        Member& m = members[i];
        const kerr::Vec4& F = force[i];
        if (F[1] * F[1] + F[2] * F[2] + F[3] * F[3] < 1e-28 || !std::isfinite(m.s.r)) continue;
        const kerr::State& s = m.s;
        const double sn = std::sin(s.theta), cs = std::cos(s.theta);
        const kerr::Vec4 u = kerr::fourVelocity(a, m.c, s);
        kerr::Vec4 low = kerr::lower(kerr::metricTensor(a, s.r, s.theta), u);
        // The force lies in the body's rest space: du_mu / dtau = F_mu over its proper
        // time, dt / u^t.
        const double dtau = dt / u[0];
        for (int k = 0; k < 4; ++k) low[k] += F[k] * dtau;

        // Back onto the unit hyperboloid: solve g^{mu nu} u_mu u_nu = -1 for u_t.
        const double sigma = kerr::metricAt(a, s.r, s.theta).sigma;
        const double delta = s.r * s.r - 2.0 * s.r + a * a;
        const double f = 2.0 * s.r / sigma;
        const double A = -(1.0 + f), B = 2.0 * f * low[1];
        const double C = delta / sigma * low[1] * low[1] + 2.0 * a / sigma * low[1] * low[3] + low[2] * low[2] / sigma +
                         low[3] * low[3] / (sigma * std::max(sn * sn, 1e-18)) + 1.0;
        const double disc = std::max(B * B - 4.0 * A * C, 0.0);
        // Of the two roots, the future-directed one has u^t = -(1 + f) u_t + f u_r > 0.
        const double r1 = (-B + std::sqrt(disc)) / (2.0 * A), r2 = (-B - std::sqrt(disc)) / (2.0 * A);
        low[0] = (-(1.0 + f) * r1 + f * low[1] > 0.0) ? r1 : r2;

        m.c.E = -low[0];
        m.c.L = low[3];
        m.s.vtheta = low[2];
        m.s.vr = 2.0 * s.r * low[0] + delta * low[1] + a * low[3];  // Sigma u^r
        m.c.Q = low[2] * low[2] + cs * cs * (a * a * (1.0 - m.c.E * m.c.E) + m.c.L * m.c.L / std::max(sn * sn, 1e-18));
    }
}

void addVelocities(const double v[3], const double u[3], double out[3]) {
    const double v2 = std::min(v[0] * v[0] + v[1] * v[1] + v[2] * v[2], 1.0 - 1e-12);
    const double gamma = 1.0 / std::sqrt(1.0 - v2);
    const double uv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    for (int k = 0; k < 3; ++k) {
        const double parallel = v2 > 0.0 ? uv / v2 * v[k] : 0.0;
        out[k] = (v[k] + parallel + (u[k] - parallel) / gamma) / (1.0 + uv);
    }
}

double stepLimit(double a, const std::vector<Member>& members, double maxStep) {
    double h = maxStep;
    for (size_t i = 0; i < members.size(); ++i) {
        const auto xi = kerr::toCartesian(a, members[i].s.r, members[i].s.theta, members[i].s.phi);
        for (size_t j = i + 1; j < members.size(); ++j) {
            const double m = members[i].mass + members[j].mass;
            if (m < kNegligibleMass) continue;
            const auto xj = kerr::toCartesian(a, members[j].s.r, members[j].s.theta, members[j].s.phi);
            const double d2 = (xj[0] - xi[0]) * (xj[0] - xi[0]) + (xj[1] - xi[1]) * (xj[1] - xi[1]) +
                              (xj[2] - xi[2]) * (xj[2] - xi[2]);
            const double eps = 0.5 * (members[i].softening + members[j].softening);
            const double d = std::sqrt(d2 + eps * eps);
            // 2 pi sqrt(d^3 / m) is the period of a circular orbit at separation d.
            h = std::min(h, 0.01 * 2.0 * kerr::kPi * std::sqrt(d * d * d / m));
        }
    }
    return std::max(h, 1e-3);
}

bool interacting(const std::vector<Member>& members) {
    int heavy = 0;
    for (const Member& m : members) heavy += m.mass >= kNegligibleMass ? 1 : 0;
    return heavy > 0 && members.size() > 1;
}

std::vector<trajectory::Path> predict(double a, std::vector<Member> members, double duration, double minSpacing,
                                      size_t maxPoints, double surface, double escapeRadius, double maxStep) {
    // Black holes: stop just outside the horizon. Stars and planets: at the surface.
    const double rh = surface > 0.0 ? surface - 0.02 : kerr::horizonRadius(a);
    const trajectory::Fate fallsIn = surface > 0.0 ? trajectory::Fate::HitsSurface : trajectory::Fate::Plunges;
    std::vector<trajectory::Path> paths(members.size());
    std::vector<bool> done(members.size(), false);
    std::vector<std::array<double, 3>> last(members.size());
    std::vector<Member> active;
    std::vector<size_t> index;
    for (size_t i = 0; i < members.size(); ++i) {
        last[i] = kerr::toCartesian(a, members[i].s.r, members[i].s.theta, members[i].s.phi);
        paths[i].points.push_back(last[i]);
        paths[i].rMin = paths[i].rMax = members[i].s.r;
    }

    double elapsed = 0.0;
    for (int guard = 0; elapsed < duration && guard < 20000; ++guard) {
        // Only members still in play pull on each other.
        active.clear();
        index.clear();
        for (size_t i = 0; i < members.size(); ++i) {
            if (!done[i]) {
                active.push_back(members[i]);
                index.push_back(i);
            }
        }
        if (active.empty()) break;
        const double h = std::min(stepLimit(a, active, maxStep), duration - elapsed);
        kick(a, active, 0.5 * h);
        for (size_t k = 0; k < active.size(); ++k) {
            Member& m = active[k];
            const double target = m.s.t + h;
            for (int g = 0; m.s.t < target - 1e-12 && g < 10000; ++g) {
                const double step = std::min(target - m.s.t, kerr::suggestedTimeStep(a, m.c, m.s, 0.01));
                kerr::stepCoordinateTime(a, m.c, m.s, step);
                kerr::normalizeAngles(m.s);
                if (!std::isfinite(m.s.r) || m.s.r < rh + 0.02 || m.s.r > escapeRadius) break;
            }
        }
        kick(a, active, 0.5 * h);
        elapsed += h;
        for (size_t k = 0; k < active.size(); ++k) {
            const size_t i = index[k];
            members[i] = active[k];
            trajectory::Path& path = paths[i];
            const kerr::State& s = members[i].s;
            path.duration = elapsed;
            if (!std::isfinite(s.r) || s.r < rh + 0.02) {
                path.fate = fallsIn;
                done[i] = true;
                continue;
            }
            if (s.r > escapeRadius) {
                path.fate = trajectory::Fate::Escapes;
                done[i] = true;
                continue;
            }
            path.rMin = std::min(path.rMin, s.r);
            path.rMax = std::max(path.rMax, s.r);
            const auto p = kerr::toCartesian(a, s.r, s.theta, s.phi);
            const double d2 = (p[0] - last[i][0]) * (p[0] - last[i][0]) + (p[1] - last[i][1]) * (p[1] - last[i][1]) +
                              (p[2] - last[i][2]) * (p[2] - last[i][2]);
            if (d2 >= minSpacing * minSpacing && path.points.size() < maxPoints) {
                path.points.push_back(p);
                last[i] = p;
            }
        }
    }
    for (size_t i = 0; i < members.size(); ++i) {
        const kerr::State& s = members[i].s;
        if (std::isfinite(s.r)) paths[i].points.push_back(kerr::toCartesian(a, s.r, s.theta, s.phi));
    }
    return paths;
}

}  // namespace gravity
