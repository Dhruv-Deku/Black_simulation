#include "eht.h"

#include "kerr.h"

#include <algorithm>
#include <cmath>
#include <thread>

namespace eht {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kStartRadius = 1000.0;  // rays start this far out (M)
constexpr double kLightRadius = 60.0;    // the flow's light is gathered inside this radius (M)
constexpr double kMaxDepth = 30.0;       // beyond this optical depth nothing behind shows

// Boyer-Lindquist metric pieces at (r, theta).
struct Metric {
    double gtt, gtp, gpp, grr;
};
Metric metric(double a, double r, double theta) {
    const double s2 = std::sin(theta) * std::sin(theta), c2 = 1.0 - s2;
    const double sigma = r * r + a * a * c2, delta = r * r - 2.0 * r + a * a;
    const double A = (r * r + a * a) * (r * r + a * a) - a * a * delta * s2;
    return {-(1.0 - 2.0 * r / sigma), -2.0 * a * r * s2 / sigma, A * s2 / sigma, sigma / delta};
}

}  // namespace

// --- Targets ----------------------------------------------------------------------------------

const std::vector<Target>& targets() {
    static const std::vector<Target> list = {
        // The spin axis points away from us (the jet's counterpart is the one we see), so
        // the inclination is 180 - 17 degrees (EHT Paper V).
        {"M87*", 6.5e9, 16.8e6, 163.0, 0.5, 42.0, 3.0, "EHT Collaboration 2019, ApJL 875, L1 and L6",
         "https://eventhorizontelescope.org/press-release-april-10-2019-astronomers-capture-first-image-black-hole"},
        // A fairly face-on view is favoured (EHT 2022, Paper V).
        {"Sgr A*", 4.3e6, 8.28e3, 30.0, 2.4, 51.8, 2.3, "EHT Collaboration 2022, ApJL 930, L12",
         "https://eventhorizontelescope.org/blog/astronomers-reveal-first-image-black-hole-heart-our-galaxy"},
    };
    return list;
}

double angularScale(double massSolar, double distancePc) {
    return kSolarLengthCm * massSolar / (distancePc * kParsecCm) / kMicroarcsecond;
}

// --- The shadow --------------------------------------------------------------------------------

std::vector<Point> criticalCurve(double a, double inclination, int samples) {
    std::vector<Point> curve;
    const double sinI = std::sin(inclination), cosI = std::cos(inclination);
    if (std::abs(a) < 1e-6) {
        for (int k = 0; k < samples; ++k) {
            const double psi = 2.0 * kPi * k / samples;
            curve.push_back({std::sqrt(27.0) * std::cos(psi), std::sqrt(27.0) * std::sin(psi)});
        }
        return curve;
    }
    // Spherical photon orbits run from the prograde to the retrograde circular one.
    const double r1 = kerr::photonOrbitRadius(a, true), r2 = kerr::photonOrbitRadius(a, false);
    std::vector<Point> upper;
    const int n = std::max(samples, 16) * 8;
    for (int k = 0; k <= n; ++k) {
        // Denser near the ends, where the curve turns round.
        const double x = 0.5 - 0.5 * std::cos(kPi * k / n);
        const double r = r1 + (r2 - r1) * x;
        const double xi = (r * r * (3.0 - r) - a * a * (r + 1.0)) / (a * (r - 1.0));
        const double eta = r * r * r * (4.0 * a * a - r * (r - 3.0) * (r - 3.0)) / (a * a * (r - 1.0) * (r - 1.0));
        const double b2 = eta + a * a * cosI * cosI - xi * xi * cosI * cosI / (sinI * sinI);
        if (b2 < 0.0) continue;
        upper.push_back({-xi / sinI, std::sqrt(b2)});
    }
    // Round: along the top, then back along the bottom.
    curve = upper;
    for (size_t k = upper.size(); k-- > 0;) curve.push_back({upper[k].alpha, -upper[k].beta});
    return curve;
}

double criticalRadius(const std::vector<Point>& curve, double psi) {
    const double dx = std::cos(psi), dy = std::sin(psi);
    double best = 0.0;
    for (size_t k = 0; k < curve.size(); ++k) {
        const Point& p = curve[k];
        const Point& q = curve[(k + 1) % curve.size()];
        // p + t (q - p) = s d, for t in [0, 1] and s > 0.
        const double ex = q.alpha - p.alpha, ey = q.beta - p.beta;
        const double det = ex * dy - ey * dx;
        if (std::abs(det) < 1e-15) continue;
        const double t = (p.beta * dx - p.alpha * dy) / det;
        const double s = (ex * p.beta - ey * p.alpha) / det;
        if (t >= 0.0 && t <= 1.0 && s > 0.0) best = std::max(best, s);
    }
    return best;
}

// --- Synchrotron radiation -------------------------------------------------------------------

namespace {

double besselI0(double x) {
    const double ax = std::abs(x);
    if (ax < 3.75) {
        const double y = (x / 3.75) * (x / 3.75);
        return 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    }
    const double y = 3.75 / ax;
    return (std::exp(ax) / std::sqrt(ax)) *
           (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 +
            y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
}

double besselI1(double x) {
    const double ax = std::abs(x);
    double result;
    if (ax < 3.75) {
        const double y = (x / 3.75) * (x / 3.75);
        result = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
    } else {
        const double y = 3.75 / ax;
        double tail = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
        tail = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * tail))));
        result = tail * std::exp(ax) / std::sqrt(ax);
    }
    return x < 0.0 ? -result : result;
}

// K0 and K1 times exp(x) (for x > 2 the exponential is factored out, so that the
// emissivity can be formed without underflow at low temperatures).
double scaledK0(double x) {
    if (x <= 2.0) {
        const double y = x * x / 4.0;
        return std::exp(x) * ((-std::log(x / 2.0) * besselI0(x)) +
               (-0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5)))))));
    }
    const double y = 2.0 / x;
    return (1.0 / std::sqrt(x)) * (1.25331414 + y * (-0.7832358e-1 + y * (0.2189568e-1 + y * (-0.1062446e-1 +
           y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))));
}

double scaledK1(double x) {
    if (x <= 2.0) {
        const double y = x * x / 4.0;
        return std::exp(x) * ((std::log(x / 2.0) * besselI1(x)) + (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579 +
               y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4))))))));
    }
    const double y = 2.0 / x;
    return (1.0 / std::sqrt(x)) * (1.25331414 + y * (0.23498619 + y * (-0.3655620e-1 + y * (0.1504268e-1 +
           y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
}

double scaledK2(double x) {
    return scaledK0(x) + 2.0 / x * scaledK1(x);
}

}  // namespace

double besselK0(double x) { return scaledK0(x) * std::exp(-x); }
double besselK1(double x) { return scaledK1(x) * std::exp(-x); }
double besselK2(double x) { return scaledK2(x) * std::exp(-x); }

double synchrotronEmissivity(double nu, double ne, double thetaE, double B, double sinThetaB) {
    if (ne <= 0.0 || thetaE < 0.05 || B <= 0.0 || sinThetaB <= 0.0) return 0.0;
    const double nuC = kElectronCharge * B / (2.0 * kPi * kElectronMass * kLight);
    const double nuS = (2.0 / 9.0) * nuC * thetaE * thetaE * sinThetaB;
    const double X = nu / nuS;
    const double cube = std::cbrt(X);
    const double shape = std::sqrt(X) + std::pow(2.0, 11.0 / 12.0) * std::pow(X, 1.0 / 6.0);
    // 1 / K2(1/thetaE) = exp(1/thetaE) / scaledK2: combined with exp(-X^1/3) in one exponent.
    const double inv = 1.0 / thetaE;
    const double exponent = inv - cube;
    if (exponent < -700.0) return 0.0;
    return ne * std::sqrt(2.0) * kPi * kElectronCharge * kElectronCharge * nuS / (3.0 * scaledK2(inv) * kLight) *
           shape * shape * std::exp(exponent);
}

double planck(double nu, double T) {
    const double x = kPlanck * nu / (kBoltzmann * T);
    return 2.0 * kPlanck * nu * nu * nu / (kLight * kLight) / std::expm1(x);
}

double brightnessTemperature(double intensity, double nu) {
    return kLight * kLight * intensity / (2.0 * kBoltzmann * nu * nu);
}

// --- The flow ---------------------------------------------------------------------------------

namespace {

// |u^r| of the geodesic plunging from the ISCO (equatorial), at r inside it.
double plungeSpeed(double a, double r) {
    const double isco = kerr::iscoRadius(a, true);
    const kerr::Constants c = kerr::circularOrbit(a, isco, true);
    const double delta = r * r - 2.0 * r + a * a, x = c.L - a * c.E;
    const double P = c.E * (r * r + a * a) - a * c.L;
    return std::sqrt(std::max(P * P - delta * (r * r + x * x), 0.0)) / (r * r);
}

}  // namespace

Plasma plasmaAt(const Flow& flow, double r, double theta, double a) {
    const double z = r * std::cos(theta), rho = std::max(r * std::sin(theta), 1e-6);
    Plasma p;
    const double h = std::max(flow.height, 0.05) * rho;
    double radial = std::pow(r, -1.1);
    const double isco = kerr::iscoRadius(a, true);
    if (r < isco) {
        // Inside the ISCO the gas falls freely: by continuity n ~ 1 / (r^2 u^r), matched to the
        // power law a quarter M inside the ISCO (at the ISCO itself the plunge starts from rest).
        const double start = std::max(isco - 0.25, 0.5 * (isco + kerr::horizonRadius(a)));
        if (r < start) {
            const double here = plungeSpeed(a, r), there = plungeSpeed(a, start);
            if (here > 0.0) radial = std::min(radial, std::pow(start, -1.1) * (start * start * there) / (r * r * here));
        }
    }
    p.ne = flow.n0 * radial * std::exp(-z * z / (2.0 * h * h));
    p.Te = flow.T0 * std::pow(r, -0.84);
    p.B = std::sqrt(8.0 * kPi * p.ne * kProtonMass * kLight * kLight / (6.0 * flow.beta * r));
    return p;
}

void flowVelocity(const Flow& flow, double a, double r, double theta, double u[4]) {
    const Metric g = metric(a, r, theta);
    double omega = 0.0, ur = 0.0;
    const double isco = kerr::iscoRadius(a, true);
    if (r >= isco) {
        omega = flow.keplerian / (std::pow(r, 1.5) + a);
    } else {
        // Plunging from the ISCO with its energy and angular momentum (equatorial geodesic).
        const kerr::Constants c = kerr::circularOrbit(a, isco, true);
        const double delta = r * r - 2.0 * r + a * a, x = c.L - a * c.E;
        const double P = c.E * (r * r + a * a) - a * c.L;
        const double R = std::max(P * P - delta * (r * r + x * x), 0.0);
        const double ut = ((r * r + a * a) * P / delta + a * x) / (r * r);
        const double up = (a * P / delta + x) / (r * r);
        omega = up / ut;
        ur = -std::sqrt(R) / (r * r);
    }
    // u^t from u.u = -1 at this (r, theta), keeping omega and u^r.
    double denom = g.gtt + 2.0 * omega * g.gtp + omega * omega * g.gpp;
    double num = -1.0 - g.grr * ur * ur;
    if (!(denom < 0.0) || !(num < 0.0)) {
        // Not possible here (deep inside, off the plane): the local ZAMO's motion.
        omega = -g.gtp / g.gpp;
        ur = 0.0;
        denom = g.gtt + 2.0 * omega * g.gtp + omega * omega * g.gpp;
        num = -1.0;
    }
    const double ut = std::sqrt(num / denom);
    u[0] = ut;
    u[1] = ur;
    u[2] = 0.0;
    u[3] = omega * ut;
}

// --- The image ----------------------------------------------------------------------------------

Ray traceRay(const View& camera, const Flow& flow, double alpha, double beta, bool withLight) {
    const double a = camera.a;
    const double sinI = std::sin(camera.inclination), cosI = std::cos(camera.inclination);
    // Bardeen: L/E = -alpha sin i, Q/E^2 = beta^2 + (alpha^2 - a^2) cos^2 i, with E = 1.
    kerr::Constants c{1.0, -alpha * sinI, beta * beta + (alpha * alpha - a * a) * cosI * cosI, 0.0};
    kerr::State s{};
    s.r = kStartRadius;
    s.theta = std::clamp(camera.inclination, 1e-4, kPi - 1e-4);
    s.vr = -std::sqrt(std::max(kerr::radialPotential(a, c, s.r), 0.0));
    // Traced back, a ray that arrives above the center (beta > 0) heads for smaller theta.
    s.vtheta = (beta > 0.0 ? -1.0 : 1.0) * std::sqrt(std::max(kerr::polarPotential(a, c, s.theta), 0.0));

    Ray ray;
    ray.minRadius = s.r;
    const double horizon = kerr::horizonRadius(a);
    const double lengthCm = kSolarLengthCm * camera.massSolar;
    double depth = 0.0;
    for (int step = 0; step < 40000; ++step) {
        const double tolerance = s.r < kLightRadius ? 0.02 : 0.05;
        const double h = kerr::photonStepSize(a, c, s, tolerance);
        const double sigmaBefore = s.r * s.r + a * a * std::cos(s.theta) * std::cos(s.theta);
        kerr::stepMino(a, c, s, h);
        kerr::enforceConstraints(a, c, s);
        kerr::normalizeAngles(s);
        ray.minRadius = std::min(ray.minRadius, s.r);
        if (!std::isfinite(s.r) || s.r < horizon * (1.0 + 1e-4)) {
            ray.captured = true;
            break;
        }
        if (s.r > kStartRadius * 1.01 && s.vr > 0.0) break;  // back out to the observer's distance
        if (!withLight || s.r > kLightRadius) continue;

        // The gas here, and the frequency it must emit at to be seen at camera.frequency.
        const double sigma = s.r * s.r + a * a * std::cos(s.theta) * std::cos(s.theta);
        const double dLambda = 0.5 * (sigma + sigmaBefore) * h * lengthCm;  // affine length, cm
        const Plasma plasma = plasmaAt(flow, s.r, s.theta, a);
        double u[4];
        flowVelocity(flow, a, s.r, s.theta, u);
        const double delta = s.r * s.r - 2.0 * s.r + a * a;
        // The photon's covariant momentum, travelling toward the camera (opposite to the
        // trace): p_t = -1, p_phi = L, p_r = -vr / Delta, p_theta = -vtheta.
        const double pr = -s.vr / delta, pth = -s.vtheta;
        const double minusPu = u[0] - c.L * u[3] - pr * u[1] - pth * u[2];
        if (!(minusPu > 0.0)) continue;
        const double g = 1.0 / minusPu;
        const double nuF = camera.frequency / g;
        const double thetaE = kBoltzmann * plasma.Te / (kElectronMass * kLight * kLight);
        const double j = synchrotronEmissivity(nuF, plasma.ne, thetaE, plasma.B, flow.sinThetaB);
        if (j <= 0.0) continue;
        const double absorption = j / planck(nuF, plasma.Te);
        const double dTau = absorption * dLambda / g;
        // Through the step, exactly for a constant source function S = j / alpha:
        // g^3 S (1 - exp(-dTau)), dimmed by what lies in front.
        const double added = dTau < 1e-6 ? g * g * j * dLambda : g * g * g * (j / absorption) * -std::expm1(-dTau);
        ray.intensity += std::exp(-depth) * added;
        depth += dTau;
        if (depth > kMaxDepth) break;
    }
    return ray;
}

Image render(const View& camera, const Flow& flow, int threads, const std::atomic<bool>* cancel) {
    Image image;
    image.size = std::max(camera.size, 1);
    image.halfWidth = camera.halfWidth;
    const size_t n = static_cast<size_t>(image.size) * image.size;
    image.intensity.assign(n, 0.0);
    image.shadow.assign(n, 0);
    std::atomic<int> nextRow{0};
    auto work = [&]() {
        for (int row = nextRow++; row < image.size; row = nextRow++) {
            if (cancel && cancel->load()) return;
            const double beta = image.halfWidth - (row + 0.5) * image.pixel();
            for (int col = 0; col < image.size; ++col) {
                const double alpha = -image.halfWidth + (col + 0.5) * image.pixel();
                const Ray ray = traceRay(camera, flow, alpha, beta);
                const size_t k = static_cast<size_t>(row) * image.size + col;
                image.intensity[k] = ray.intensity;
                image.shadow[k] = ray.captured ? 1 : 0;
            }
        }
    };
    std::vector<std::thread> pool;
    for (int t = 1; t < std::max(threads, 1); ++t) pool.emplace_back(work);
    work();
    for (std::thread& t : pool) t.join();
    return image;
}

double flux(const Image& image, double thetaGUas) {
    const double solidAngle = std::pow(image.pixel() * thetaGUas * kMicroarcsecond, 2.0);
    double sum = 0.0;
    for (double i : image.intensity) sum += i;
    return sum * solidAngle / kJansky;
}

Image blur(const Image& image, double fwhm) {
    Image out = image;
    const double sigma = fwhm / 2.354820045 / image.pixel();
    if (sigma < 0.3) return out;
    const int reach = static_cast<int>(std::ceil(3.5 * sigma));
    std::vector<double> kernel(2 * reach + 1);
    double norm = 0.0;
    for (int k = -reach; k <= reach; ++k) norm += kernel[k + reach] = std::exp(-0.5 * k * k / (sigma * sigma));
    for (double& w : kernel) w /= norm;
    const int n = image.size;
    std::vector<double> tmp(image.intensity.size(), 0.0);
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            double sum = 0.0;
            for (int k = -reach; k <= reach; ++k) {
                const int xx = x + k;
                if (xx >= 0 && xx < n) sum += kernel[k + reach] * image.intensity[static_cast<size_t>(y) * n + xx];
            }
            tmp[static_cast<size_t>(y) * n + x] = sum;
        }
    }
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            double sum = 0.0;
            for (int k = -reach; k <= reach; ++k) {
                const int yy = y + k;
                if (yy >= 0 && yy < n) sum += kernel[k + reach] * tmp[static_cast<size_t>(yy) * n + x];
            }
            out.intensity[static_cast<size_t>(y) * n + x] = sum;
        }
    }
    return out;
}

Point curveCenter(const std::vector<Point>& curve) {
    Point c{0.0, 0.0};
    for (const Point& p : curve) {
        c.alpha += p.alpha / curve.size();
        c.beta += p.beta / curve.size();
    }
    return c;
}

double ringDiameter(const Image& image, Point center, double radius) {
    const int n = image.size;
    // Pixels: x right from the left edge, y down from the top.
    const double cx = (center.alpha + image.halfWidth) / image.pixel();
    const double cy = (image.halfWidth - center.beta) / image.pixel();
    auto sample = [&](double x, double y) {
        x -= 0.5;
        y -= 0.5;
        const int x0 = static_cast<int>(std::floor(x)), y0 = static_cast<int>(std::floor(y));
        if (x0 < 0 || y0 < 0 || x0 + 1 >= n || y0 + 1 >= n) return 0.0;
        const double fx = x - x0, fy = y - y0;
        auto at = [&](int i, int j) { return image.intensity[static_cast<size_t>(j) * n + i]; };
        return (1 - fx) * (1 - fy) * at(x0, y0) + fx * (1 - fy) * at(x0 + 1, y0) + (1 - fx) * fy * at(x0, y0 + 1) +
               fx * fy * at(x0 + 1, y0 + 1);
    };
    const double from = 0.4 * radius / image.pixel(), to = 2.0 * radius / image.pixel();
    double total = 0.0;
    const int directions = 72;
    for (int k = 0; k < directions; ++k) {
        const double psi = 2.0 * kPi * k / directions;
        double bestR = from, best = -1.0;
        for (double r = from; r < to; r += 0.05) {
            const double v = sample(cx + r * std::cos(psi), cy - r * std::sin(psi));
            if (v > best) {
                best = v;
                bestR = r;
            }
        }
        total += bestR;
    }
    return 2.0 * total / directions * image.pixel();
}

double calibrate(const View& camera, Flow flow, double targetJy, double thetaGUas, int threads,
                 const std::atomic<bool>* cancel) {
    View small = camera;
    small.size = 48;
    auto fluxAt = [&](double n0) {
        flow.n0 = n0;
        return flux(render(small, flow, threads, cancel), thetaGUas);
    };
    double n0 = flow.n0, f0 = fluxAt(n0);
    // Find a density with some light first.
    for (int k = 0; k < 8 && f0 <= 0.0; ++k) f0 = fluxAt(n0 *= 10.0);
    if (f0 <= 0.0) return n0;
    double n1 = n0 * targetJy / f0, f1 = fluxAt(n1);
    for (int k = 0; k < 8; ++k) {
        if (cancel && cancel->load()) break;
        if (f1 <= 0.0 || std::abs(f1 / targetJy - 1.0) < 0.01) break;
        // Secant in log n0 against log flux.
        const double slope = (std::log(f1) - std::log(f0)) / (std::log(n1) - std::log(n0));
        const double next = std::exp(std::log(n1) + (std::log(targetJy) - std::log(f1)) / (std::abs(slope) > 0.05 ? slope : 1.0));
        n0 = n1;
        f0 = f1;
        n1 = std::clamp(next, n0 / 30.0, n0 * 30.0);
        f1 = fluxAt(n1);
    }
    return n1;
}

void afmhot(double x, float rgb[3]) {
    x = std::clamp(x, 0.0, 1.0);
    rgb[0] = static_cast<float>(std::clamp(2.0 * x, 0.0, 1.0));
    rgb[1] = static_cast<float>(std::clamp(2.0 * x - 0.5, 0.0, 1.0));
    rgb[2] = static_cast<float>(std::clamp(2.0 * x - 1.0, 0.0, 1.0));
}

}  // namespace eht
