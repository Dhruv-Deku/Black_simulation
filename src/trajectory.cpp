#include "trajectory.h"

#include "central.h"

#include <algorithm>
#include <cmath>

namespace trajectory {

namespace {

double dot(const Point& a, const Point& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double length(const Point& a) {
    return std::sqrt(dot(a, a));
}

double angleBetween(const Point& a, const Point& b) {
    double c = dot(a, b) / std::max(length(a) * length(b), 1e-300);
    return std::acos(std::clamp(c, -1.0, 1.0));
}

Point position(double a, const kerr::State& s) {
    return kerr::toCartesian(a, s.r, s.theta, s.phi);
}

// First zero of R between r0 and r1 (either order), searched on a log grid from r0,
// then refined by bisection. Returns false if R stays non-negative.
bool firstRoot(double a, const kerr::Constants& c, double r0, double r1, double& root) {
    const int n = 4000;
    const double ratio = std::pow(r1 / r0, 1.0 / n);
    double prev = r0;
    for (int i = 1; i <= n; ++i) {
        double r = r0 * std::pow(ratio, i);
        if (kerr::radialPotential(a, c, r) < 0.0) {
            double lo = prev, hi = r;  // R(lo) >= 0 > R(hi)
            for (int k = 0; k < 80; ++k) {
                double mid = 0.5 * (lo + hi);
                (kerr::radialPotential(a, c, mid) >= 0.0 ? lo : hi) = mid;
            }
            root = 0.5 * (lo + hi);
            return true;
        }
        prev = r;
    }
    return false;
}

}  // namespace

const char* fateName(Fate fate) {
    switch (fate) {
        case Fate::Bound: return "Bound orbit";
        case Fate::Circular: return "Circular orbit";
        case Fate::Plunges: return "Falls into the black hole";
        case Fate::Escapes: return "Escapes to infinity";
        case Fate::HitsSurface: return "Hits the surface";
    }
    return "?";
}

const char* rayEndName(RayEnd end) {
    switch (end) {
        case RayEnd::Sky: return "the distant sky";
        case RayEnd::Horizon: return "the black hole (the shadow)";
        case RayEnd::Disk: return "the accretion disk";
        case RayEnd::Body: return "an orbiting body";
        case RayEnd::Unfinished: return "nowhere yet (ran out of steps)";
        case RayEnd::Surface: return "the surface of the central object";
        case RayEnd::Lens: return "a small black hole (captured)";
    }
    return "?";
}

double effectivePotential(double a, double L, double Q, double r) {
    // R(r) = alpha E^2 - 2 beta E + gamma, with
    // alpha = (r^2 + a^2)^2 - a^2 Delta, beta = 2 a L r, gamma = a^2 L^2 - Delta (r^2 + L^2 + Q).
    const double delta = r * r - 2.0 * r + a * a;
    const double rra = r * r + a * a;
    const double alpha = rra * rra - a * a * delta;
    const double beta = 2.0 * a * L * r;
    const double gamma = a * a * L * L - delta * (r * r + L * L + Q);
    const double disc = std::max(beta * beta - alpha * gamma, 0.0);
    return (beta + std::sqrt(disc)) / alpha;
}

Analysis classify(double a, const kerr::Constants& c, const kerr::State& s, double surface) {
    Analysis out;
    out.inclinationDeg = std::atan2(std::sqrt(std::max(c.Q, 0.0)), std::abs(c.L)) * 180.0 / kerr::kPi;
    // The innermost place an orbit can reach: the horizon, or the surface.
    const double rh = central::boundaryRadius(a, surface);
    const Fate fallsIn = surface > 0.0 ? Fate::HitsSurface : Fate::Plunges;
    const double r0 = s.r;

    // Nudge off the starting radius so a body exactly at a turning point (R = 0)
    // does not count it as the next one.
    const double eps = 1e-7 * r0;
    out.hasInner = r0 - eps > rh && firstRoot(a, c, r0 - eps, rh * (1.0 + 1e-9), out.rInner);
    out.hasOuter = firstRoot(a, c, r0 + eps, std::max(1e5, 100.0 * r0), out.rOuter);

    // Which way it is heading; at a turning point, the way the potential pushes it.
    const double h = 1e-6 * r0;
    const double slope = (kerr::radialPotential(a, c, r0 + h) - kerr::radialPotential(a, c, r0 - h)) / (2 * h);
    const double scale = r0 * r0 * r0 * r0;
    const bool atTurn = s.vr * s.vr < 1e-12 * scale;
    const bool outward = atTurn ? slope > 0.0 : s.vr > 0.0;

    if (atTurn && std::abs(slope) < 1e-7 * scale / r0) {
        out.fate = Fate::Circular;
        out.stableCircular = out.hasInner && out.hasOuter && out.rOuter - out.rInner < 1e-4 * r0;
        return out;
    }
    if (out.hasInner && out.hasOuter) out.fate = Fate::Bound;
    else if (!outward) out.fate = out.hasInner ? Fate::Escapes : fallsIn;
    else out.fate = out.hasOuter ? fallsIn : Fate::Escapes;
    return out;
}

Path predict(double a, kerr::Constants c, kerr::State s, double maxTime, double minSpacing, size_t maxPoints,
             double surface, double escapeRadius) {
    Path path;
    // Black holes: stop just outside the horizon. Stars and planets: at the surface.
    const double rh = surface > 0.0 ? surface - 0.02 : kerr::horizonRadius(a);
    const double t0 = s.t;
    path.rMin = path.rMax = s.r;
    Point last = position(a, s);
    path.points.push_back(last);

    Point prevPos = last;
    double swept = 0.0;  // angle swept around the hole since the last periapsis
    double prevVr = s.vr;
    double firstPeriapsisSweep = -1.0;

    for (int guard = 0; guard < 400000 && s.t - t0 < maxTime; ++guard) {
        double h = std::min(kerr::suggestedTimeStep(a, c, s, 0.01), maxTime - (s.t - t0));
        kerr::stepCoordinateTime(a, c, s, h);
        kerr::normalizeAngles(s);
        if (!std::isfinite(s.r) || s.r < rh + 0.02) {
            path.fate = surface > 0.0 ? Fate::HitsSurface : Fate::Plunges;
            break;
        }
        if (s.r > escapeRadius) {
            path.fate = Fate::Escapes;
            break;
        }
        path.rMin = std::min(path.rMin, s.r);
        path.rMax = std::max(path.rMax, s.r);

        Point p = position(a, s);
        swept += angleBetween(prevPos, p);
        prevPos = p;
        // Periapsis: the radial velocity turns from inward to outward.
        if (prevVr < 0.0 && s.vr >= 0.0) {
            if (++path.periapses >= 2 && firstPeriapsisSweep < 0.0) firstPeriapsisSweep = swept;
            swept = 0.0;
        }
        prevVr = s.vr;

        if (std::sqrt((p[0] - last[0]) * (p[0] - last[0]) + (p[1] - last[1]) * (p[1] - last[1]) +
                      (p[2] - last[2]) * (p[2] - last[2])) >= minSpacing &&
            path.points.size() < maxPoints) {
            path.points.push_back(p);
            last = p;
        }
    }
    path.points.push_back(position(a, s));
    path.duration = s.t - t0;
    if (firstPeriapsisSweep > 0.0) path.precessionDeg = (firstPeriapsisSweep - 2.0 * kerr::kPi) * 180.0 / kerr::kPi;
    return path;
}

bool lensKick(double a, const Point& center, double mass, bool backward, kerr::Constants& c, kerr::State& s) {
    const double r = s.r, th = s.theta;
    const Point x = position(a, s);
    const kerr::Metric m = kerr::metricAt(a, r, th);
    const double sn = std::max(std::abs(std::sin(th)), 1e-9), cs = std::cos(th);
    const double alpha = std::sqrt(m.sigma * m.delta / m.A);  // ZAMO lapse
    const double varpi = std::sqrt(m.A / m.sigma) * sn;       // circumference radius
    const double sd = std::sqrt(m.sigma * m.delta), sq = std::sqrt(m.sigma);
    // The photon in the local ZAMO frame: its energy, and its direction of travel
    // forward in time (on the ZAMO's r, theta, phi axes, then in Cartesian).
    const double eLocal = (c.E - m.omega * c.L) / alpha;
    const double n[3] = {s.vr / (eLocal * sd), s.vtheta / (eLocal * sq), c.L / (eLocal * varpi)};
    const kerr::LocalAxes ax = kerr::localAxes(a, r, th, s.phi);
    Point dir{};
    for (int i = 0; i < 3; ++i) dir[i] = n[0] * ax.r[i] + n[1] * ax.theta[i] + n[2] * ax.phi[i];
    const double norm = std::max(length(dir), 1e-300);
    for (double& v : dir) v /= norm;
    // Impact parameter: the part of the separation across the direction of travel.
    const Point toLens = {center[0] - x[0], center[1] - x[1], center[2] - x[2]};
    const double along = dot(toLens, dir);
    Point perp = {toLens[0] - along * dir[0], toLens[1] - along * dir[1], toLens[2] - along * dir[2]};
    const double b = length(perp);
    if (b < 3.0 * std::sqrt(3.0) * mass) return false;
    const double bend = std::min(4.0 * mass / b, 0.5 * kerr::kPi);
    for (double& v : perp) v /= b;
    // Forward in time the ray turns toward the lens. Traced backwards we hold the later
    // direction and want the earlier one, which is turned the other way (so that the
    // backward-going ray, too, bends toward the lens).
    const double turn = backward ? -bend : bend;
    Point bent;
    for (int i = 0; i < 3; ++i) bent[i] = std::cos(turn) * dir[i] + std::sin(turn) * perp[i];
    const double nr = dot(bent, ax.r), nt = dot(bent, ax.theta), np = dot(bent, ax.phi);
    c.L = eLocal * np * varpi;
    c.E = alpha * eLocal + m.omega * c.L;
    s.vr = eLocal * nr * sd;
    s.vtheta = eLocal * nt * sq;
    c.Q = s.vtheta * s.vtheta + cs * cs * (c.L * c.L / (sn * sn) - a * a * c.E * c.E);
    return true;
}

RayReport traceRay(double a, double r, double theta, double phi, const kerr::Frame& frame, const double look[3],
                   const RayScene& scene) {
    RayReport out;
    kerr::Constants c;
    kerr::State s;
    kerr::photonFromFrame(a, r, theta, phi, frame, look, c, s);
    out.energy = c.E;
    const double invE = std::abs(c.E) > 1e-12 ? 1.0 / c.E : 0.0;
    out.gCamera = invE;
    out.b = c.L * invE;
    out.q = c.Q * invE * invE;
    out.rMin = r;

    const double rh = kerr::horizonRadius(a);
    const bool star = scene.surface > 0.0;
    const double escape = std::max(scene.escapeRadius, 1.5 * r);
    Point x = position(a, s);
    out.path.push_back(x);
    double swept = 0.0;

    // Inside the horizon every ray, followed back in time, rises to the horizon. If
    // it reaches it with E - Omega_H L < 0 (P(r+) < 0) it crossed into the hole from
    // the far side of the eternal Kerr spacetime, through which nothing is seen.
    if (!star && r < rh && 2.0 * rh * c.E - a * c.L < 0.0) {
        out.end = RayEnd::Horizon;
        out.rEnd = r;
        return out;
    }

    auto lerp = [](const kerr::State& p, const kerr::State& q, double f) {
        kerr::State o;
        o.t = p.t + f * (q.t - p.t);
        o.r = p.r + f * (q.r - p.r);
        o.theta = p.theta + f * (q.theta - p.theta);
        o.phi = p.phi + f * (q.phi - p.phi);
        o.vr = p.vr + f * (q.vr - p.vr);
        o.vtheta = p.vtheta + f * (q.vtheta - p.vtheta);
        return o;
    };

    // Each small black hole bends a ray once, at its closest approach (thin lens);
    // after the kink the ray, now heading slightly toward it, must not be kicked again.
    std::vector<bool> lensed(scene.lenses.size(), false);

    // March backwards in time along the photon's path (negative Mino-time steps).
    for (int i = 0; i < 200000; ++i) {
        const kerr::State prev = s;
        kerr::stepMino(a, c, s, -kerr::photonStepSize(a, c, s, 0.01));
        kerr::normalizeAngles(s);
        kerr::enforceConstraints(a, c, s);
        Point x1 = position(a, s);

        // A small black hole passed within this step: cut the step at the closest
        // approach and bend the ray there (thin lens), or end it if it is captured.
        // (Not across a pole, where the angles jump.)
        if (!scene.lenses.empty() && std::abs(s.phi - prev.phi) < 1.0) {
            const Point d = {x1[0] - x[0], x1[1] - x[1], x1[2] - x[2]};
            const double dd = std::max(dot(d, d), 1e-300);
            double bestU = 2.0;
            int lens = -1;
            for (size_t k = 0; k < scene.lenses.size(); ++k) {
                const Point& lc = scene.lenses[k].center;
                const double u = -((x[0] - lc[0]) * d[0] + (x[1] - lc[1]) * d[1] + (x[2] - lc[2]) * d[2]) / dd;
                if (!lensed[k] && u > 0.0 && u < 1.0 && u < bestU) {
                    bestU = u;
                    lens = static_cast<int>(k);
                }
            }
            if (lens >= 0) {
                s = lerp(prev, s, bestU);
                x1 = position(a, s);
                if (!lensKick(a, scene.lenses[lens].center, scene.lenses[lens].mass, true, c, s)) {
                    out.end = RayEnd::Lens;
                    out.body = lens;
                    out.rEnd = s.r;
                    out.delay = -s.t;
                    out.path.push_back(x1);
                    break;
                }
                lensed[lens] = true;
                ++out.lensings;
            }
        }

        out.rMin = std::min(out.rMin, s.r);
        swept += angleBetween(x, x1);

        // The surface of a star or planet within this step (as a fraction of it).
        double tSurface = 2.0;
        if (star && s.r <= scene.surface) tSurface = (prev.r - scene.surface) / std::max(prev.r - s.r, 1e-300);

        // Disk crossing: theta passes through pi/2 inside the disk (strictly: a camera
        // sitting exactly in the plane is not on the disk).
        if (scene.disk && (prev.theta - 0.5 * kerr::kPi) * (s.theta - 0.5 * kerr::kPi) < 0.0) {
            double f = (prev.theta - 0.5 * kerr::kPi) / (prev.theta - s.theta);
            double rc = prev.r + f * (s.r - prev.r);
            if (f < tSurface && rc >= scene.diskInner && rc <= scene.diskOuter) {
                // Gas on a circular orbit: g = gCam / (u^t (1 - b Omega)).
                double omega = kerr::circularOmega(a, rc, scene.diskPrograde);
                double gtt = -(1.0 - 2.0 / rc), gtp = -2.0 * a / rc, gpp = rc * rc + a * a + 2.0 * a * a / rc;
                double norm = -(gtt + 2.0 * gtp * omega + gpp * omega * omega);
                if (norm > 0.0) {
                    out.end = RayEnd::Disk;
                    out.rEnd = rc;
                    out.g = std::sqrt(norm) / (c.E - c.L * omega);
                    out.delay = -(prev.t + f * (s.t - prev.t));
                    out.path.push_back(x1);
                    break;
                }
            }
        }

        // Body hit: segment against each sphere.
        int hitBody = -1;
        double bestT = 2.0;
        const Point d = {x1[0] - x[0], x1[1] - x[1], x1[2] - x[2]};
        for (size_t k = 0; k < scene.spheres.size(); ++k) {
            const Sphere& sp = scene.spheres[k];
            const Point m = {x[0] - sp.center[0], x[1] - sp.center[1], x[2] - sp.center[2]};
            double dd = std::max(dot(d, d), 1e-300), bq = dot(m, d), cq = dot(m, m) - sp.radius * sp.radius;
            if (cq > 0.0 && bq > 0.0) continue;
            double disc = bq * bq - dd * cq;
            if (disc < 0.0) continue;
            double t = cq < 0.0 ? 0.0 : (-bq - std::sqrt(disc)) / dd;
            if (t <= 1.0 && t < bestT && t < tSurface) {
                bestT = t;
                hitBody = static_cast<int>(k);
            }
        }
        if (hitBody >= 0) {
            out.end = RayEnd::Body;
            out.body = hitBody;
            out.rEnd = prev.r + bestT * (s.r - prev.r);
            out.delay = -(prev.t + bestT * (s.t - prev.t));
            out.path.push_back(x1);
            break;
        }

        if (tSurface <= 1.0) {
            // Light leaving the rotating surface: g = 1 / (u^t (E - Omega L)), and the
            // emission angle from the normal in the surface's own frame,
            // cos = (p . e_r) / (-p . u) = g |v_r| / sqrt(Sigma Delta).
            const kerr::State hit = lerp(prev, s, tSurface);
            const double ut = central::rotatingSurfaceUt(a, scene.surface, hit.theta, scene.surfaceOmega);
            const kerr::Metric m = kerr::metricAt(a, scene.surface, hit.theta);
            out.end = RayEnd::Surface;
            out.rEnd = scene.surface;
            out.g = ut > 0.0 ? 1.0 / (ut * (c.E - scene.surfaceOmega * c.L)) : 0.0;
            // |v_r| from the null condition at the surface itself (interpolating it
            // between the steps is not accurate enough).
            const double vr = std::sqrt(std::max(kerr::radialPotential(a, c, scene.surface), 0.0));
            out.emissionCos = std::min(out.g * vr / std::sqrt(m.sigma * m.delta), 1.0);
            out.surfaceTheta = hit.theta;
            out.surfacePhi = hit.phi;
            out.delay = -hit.t;
            out.path.push_back(position(a, hit));
            break;
        }

        x = x1;
        if (out.path.size() < 20000) out.path.push_back(x1);
        // Sinking back toward the horizon from outside: it came out of the past
        // horizon, the black hole's shadow.
        if (!star && s.r < rh + 0.01 && s.r < prev.r) {
            out.end = RayEnd::Horizon;
            out.rEnd = s.r;
            out.delay = -s.t;
            break;
        }
        if (s.r > escape && s.r > prev.r) {
            out.end = RayEnd::Sky;
            out.rEnd = s.r;
            out.delay = -s.t;
            const Point dir = {x1[0] - out.path[out.path.size() - 2][0], x1[1] - out.path[out.path.size() - 2][1],
                               x1[2] - out.path[out.path.size() - 2][2]};
            const double n = std::max(length(dir), 1e-300);
            out.skyDirection = {dir[0] / n, dir[1] / n, dir[2] / n};
            break;
        }
    }
    out.turns = swept / (2.0 * kerr::kPi);
    return out;
}

}  // namespace trajectory
