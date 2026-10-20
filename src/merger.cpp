#include "merger.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace merger {

// --- Inspiral -----------------------------------------------------------------------------

double chirpMass(double m1, double m2) {
    return std::pow(m1 * m2, 0.6) / std::pow(m1 + m2, 0.2);
}

double coalescenceTime(double a, double m1, double m2) {
    return 5.0 * std::pow(a, 4) / (256.0 * m1 * m2 * (m1 + m2));
}

double separationAfter(double a, double t, double m1, double m2) {
    const double left = 1.0 - t / coalescenceTime(a, m1, m2);
    return left > 0.0 ? a * std::pow(left, 0.25) : 0.0;
}

double chirpFrequency(double tau, double mc) {
    return std::pow(5.0 / (256.0 * std::max(tau, 1e-300)), 0.375) * std::pow(mc, -0.625) / kPi;
}

double gwLuminosity(double a, double m1, double m2) {
    const double M = m1 + m2, eta = m1 * m2 / (M * M);
    return 32.0 / 5.0 * eta * eta * std::pow(M / a, 5);
}

glm::dvec3 radiationReaction(const glm::dvec3& x, const glm::dvec3& v, double m1, double m2) {
    const double M = m1 + m2, eta = m1 * m2 / (M * M);
    const double r = glm::length(x);
    const glm::dvec3 n = x / r;
    const double rdot = glm::dot(n, v), v2 = glm::dot(v, v);
    return 8.0 / 5.0 * eta * M * M / (r * r * r) *
           ((3.0 * v2 + 17.0 / 3.0 * M / r) * rdot * n - (v2 + 3.0 * M / r) * v);
}

Strain strain(const std::array<double, 6>& iddot, const glm::dvec3& nIn, double distance) {
    const glm::dvec3 n = glm::normalize(nIn);
    const double theta = std::acos(std::clamp(n.z, -1.0, 1.0));
    const double phi = std::atan2(n.y, n.x);
    const glm::dvec3 p(std::cos(theta) * std::cos(phi), std::cos(theta) * std::sin(phi), -std::sin(theta));
    const glm::dvec3 q(-std::sin(phi), std::cos(phi), 0.0);
    const double h[3][3] = {{iddot[0], iddot[3], iddot[4]}, {iddot[3], iddot[1], iddot[5]}, {iddot[4], iddot[5], iddot[2]}};
    double pp = 0.0, qq = 0.0, pq = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            pp += p[i] * p[j] * h[i][j];
            qq += q[i] * q[j] * h[i][j];
            pq += p[i] * q[j] * h[i][j];
        }
    // h_ij = (2 / D) Iddot_ij; h+ = (h_pp - h_qq) / 2, hx = h_pq.
    return {(pp - qq) / distance, 2.0 * pq / distance};
}

std::array<double, 6> quadrupoleSecondDerivative(const std::vector<glm::dvec3>& x, const std::vector<glm::dvec3>& v,
                                                 const std::vector<glm::dvec3>& a, const std::vector<double>& m) {
    std::array<double, 6> out{};
    const int pairs[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
    for (size_t k = 0; k < x.size(); ++k) {
        for (int c = 0; c < 6; ++c) {
            const int i = pairs[c][0], j = pairs[c][1];
            out[c] += m[k] * (2.0 * v[k][i] * v[k][j] + x[k][i] * a[k][j] + x[k][j] * a[k][i]);
        }
    }
    return out;
}

// --- The stars -----------------------------------------------------------------------------

double polytropeK(double radius) {
    return 2.0 * radius * radius / kPi;
}

double polytropeCentralDensity(double mass, double radius) {
    return kPi * mass / (4.0 * radius * radius * radius);
}

double polytropeDensity(double mass, double radius, double r) {
    if (r >= radius) return 0.0;
    const double xi = kPi * r / radius;
    return polytropeCentralDensity(mass, radius) * (xi < 1e-8 ? 1.0 : std::sin(xi) / xi);
}

double polytropeEnclosed(double mass, double radius, double r) {
    if (r >= radius) return mass;
    const double xi = kPi * r / radius;
    return mass * (std::sin(xi) - xi * std::cos(xi)) / kPi;
}

double polytropeRadiusOfFraction(double radius, double f) {
    f = std::clamp(f, 0.0, 1.0);
    double lo = 0.0, hi = radius;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        (polytropeEnclosed(1.0, radius, mid) < f ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

std::vector<glm::dvec3> starParticles(int n, double radius) {
    std::vector<glm::dvec3> out;
    if (n <= 0) return out;
    // Face-centred cubic lattice (4 points per cube of side A), a few percent denser than
    // needed so the n points nearest the centre fill the unit sphere.
    const double A = std::cbrt(4.0 * (4.0 * kPi / 3.0) / (1.08 * n));
    const int span = static_cast<int>(std::ceil(1.15 / A)) + 1;
    const glm::dvec3 basis[4] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    std::vector<glm::dvec3> lattice;
    for (int i = -span; i <= span; ++i)
        for (int j = -span; j <= span; ++j)
            for (int k = -span; k <= span; ++k)
                for (const glm::dvec3& b : basis) {
                    // Offset by a quarter cell so no point sits exactly at the centre.
                    const glm::dvec3 p = (glm::dvec3(i, j, k) + b + 0.25) * A;
                    if (glm::dot(p, p) < 1.3) lattice.push_back(p);
                }
    std::sort(lattice.begin(), lattice.end(),
              [](const glm::dvec3& a, const glm::dvec3& b) { return glm::dot(a, a) < glm::dot(b, b); });
    lattice.resize(std::min<size_t>(lattice.size(), static_cast<size_t>(n)));
    out.reserve(lattice.size());
    const double count = static_cast<double>(lattice.size());
    for (size_t k = 0; k < lattice.size(); ++k) {
        // The k-th point holds the k-th equal share of the mass.
        const double f = (static_cast<double>(k) + 0.5) / count;
        const double len = glm::length(lattice[k]);
        out.push_back(lattice[k] / std::max(len, 1e-12) * polytropeRadiusOfFraction(radius, f));
    }
    return out;
}

double temperatureFromThermalEnergy(double uThermal) {
    // m_n c^2 / k = 939.56542 MeV / 8.617333e-11 MeV/K.
    return 2.0 / 3.0 * 1.090327e13 * std::max(uThermal, 0.0);
}

double thermalEnergyFromTemperature(double kelvin) {
    return 1.5 * kelvin / 1.090327e13;
}

// --- Collapse --------------------------------------------------------------------------------

BlackHole collapse(const std::vector<glm::dvec3>& x, const std::vector<glm::dvec3>& v, const std::vector<double>& m,
                   const std::vector<double>& phi, const std::vector<bool>& alive, const glm::dvec3& center,
                   std::vector<bool>& swallowed) {
    const size_t n = x.size();
    swallowed.assign(n, false);
    std::vector<size_t> order;
    for (size_t i = 0; i < n; ++i)
        if (alive[i]) order.push_back(i);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return glm::dot(x[a] - center, x[a] - center) < glm::dot(x[b] - center, x[b] - center);
    });
    BlackHole bh;
    if (order.empty()) return bh;

    // The core's velocity: the innermost tenth of the mass.
    double total = 0.0;
    for (size_t i : order) total += m[i];
    glm::dvec3 coreMomentum(0.0);
    double coreMass = 0.0;
    for (size_t i : order) {
        if (coreMass > 0.1 * total) break;
        coreMass += m[i];
        coreMomentum += m[i] * v[i];
    }
    const glm::dvec3 coreVelocity = coreMomentum / std::max(coreMass, 1e-300);

    double enclosed = 0.0;
    for (size_t i : order) {
        enclosed += m[i];
        const glm::dvec3 dx = x[i] - center, dv = v[i] - coreVelocity;
        const bool bound = 0.5 * glm::dot(dv, dv) + phi[i] < 0.0;
        const double j = glm::length(glm::cross(dx, dv));
        if (!bound) {
            bh.unboundMass += m[i];
        } else if (j < 2.0 * std::sqrt(3.0) * enclosed) {
            swallowed[i] = true;
            bh.mass += m[i];
            bh.position += m[i] * x[i];
            bh.velocity += m[i] * v[i];
        } else {
            bh.diskMass += m[i];
        }
    }
    if (bh.mass <= 0.0) return bh;
    bh.position /= bh.mass;
    bh.velocity /= bh.mass;
    glm::dvec3 J(0.0);
    for (size_t i = 0; i < n; ++i) {
        if (swallowed[i]) J += m[i] * glm::cross(x[i] - bh.position, v[i] - bh.velocity);
    }
    const double Jmag = glm::length(J);
    bh.spin = std::min(Jmag / (bh.mass * bh.mass), 0.998);
    bh.axis = Jmag > 0.0 ? J / Jmag : glm::dvec3(0, 0, 1);
    return bh;
}

Ringdown ringdown(double mass, double spin) {
    const double s = 1.0 - std::clamp(spin, 0.0, 0.9999);
    Ringdown r;
    r.omega = (1.5251 - 1.1568 * std::pow(s, 0.1292)) / mass;
    const double Q = 0.7 + 1.4187 * std::pow(s, -0.4990);
    r.tau = 2.0 * Q / r.omega;
    return r;
}

double horizonRadius(double mass, double spin) {
    const double a = std::clamp(spin, 0.0, 1.0) * mass;
    return mass + std::sqrt(std::max(mass * mass - a * a, 0.0));
}

// --- The debris cloud --------------------------------------------------------------------------

double debrisRadiusKm(const DebrisCloud& c, double t) {
    return c.launchKm + c.edgeSpeed * kLightCm * 1e-5 * std::max(t, 0.0);
}

double debrisTemperature(const DebrisCloud& c, double t) {
    return c.launchTemperature * c.launchKm / debrisRadiusKm(c, t);
}

double debrisOpticalDepth(const DebrisCloud& c, double t) {
    const double rCm = debrisRadiusKm(c, t) * 1e5;
    return 3.0 * c.opacity * c.mass * kSunGrams / (4.0 * kPi * rCm * rCm);
}

// --- Kilonova --------------------------------------------------------------------------------

std::array<EjectaComponent, 3> gw170817Ejecta() {
    return {{
        {"blue", 0.020, 0.256, 0.5, 674.0},
        {"purple", 0.047, 0.152, 3.0, 1308.0},
        {"red", 0.011, 0.137, 10.0, 3745.0},
    }};
}

double rProcessHeating(double t) {
    const double x = 0.5 - std::atan((t - 1.3) / 0.11) / kPi;
    return 4e18 * std::pow(std::max(x, 0.0), 1.3);
}

double thermalization(double t) {
    const double a = 0.56, b = 0.17, d = 0.74;
    const double days = std::max(t, 1e-6) / kDay;
    const double bt = 2.0 * b * std::pow(days, d);
    return 0.36 * (std::exp(-a * days) + std::log1p(bt) / bt);
}

Kilonova::Kilonova(const EjectaComponent& c) : m_c(c) {
    const double sigmaSB = 5.670374e-5;
    const double beta = 13.7;
    const double M = c.mass * kSunGrams;
    const double v = c.velocity * kLightCm;
    m_td = std::sqrt(2.0 * c.opacity * M / (beta * v * kLightCm));
    m_logT0 = -3.0;
    m_dLogT = 0.004;
    const int count = static_cast<int>((8.0 - m_logT0) / m_dLogT) + 1;
    m_L.resize(count);
    m_T.resize(count);
    m_R.resize(count);
    auto source = [&](double t) { return 2.0 * M * rProcessHeating(t) * thermalization(t) * t / (m_td * m_td); };
    double y = 0.5 * source(std::pow(10.0, m_logT0)) * std::pow(10.0, m_logT0);
    double prevT = std::pow(10.0, m_logT0);
    for (int k = 0; k < count; ++k) {
        const double t = std::pow(10.0, m_logT0 + k * m_dLogT);
        if (k > 0) {
            // y(t) = exp(-t^2/td^2) * integral of source(t') exp(t'^2/td^2): advanced stably,
            // damping what came before and adding this step (trapezoid).
            const double decay = std::exp(-(t * t - prevT * prevT) / (m_td * m_td));
            y = y * decay + 0.5 * (t - prevT) * (source(prevT) * decay + source(t));
        }
        prevT = t;
        m_L[k] = y;
        const double rEjecta = v * t;
        double T = std::pow(y / (4.0 * kPi * sigmaSB * rEjecta * rEjecta), 0.25);
        double R = rEjecta;
        if (T < c.temperatureFloor) {
            T = c.temperatureFloor;
            R = std::sqrt(y / (4.0 * kPi * sigmaSB * std::pow(T, 4)));
        }
        m_T[k] = T;
        m_R[k] = R;
        if (y > m_peakLuminosity) {
            m_peakLuminosity = y;
            m_peakTime = t;
        }
    }
}

double Kilonova::sample(const std::vector<double>& table, double t) const {
    const double u = (std::log10(std::max(t, 1e-30)) - m_logT0) / m_dLogT;
    if (u <= 0.0) return table.front();
    const int k = static_cast<int>(u);
    if (k + 1 >= static_cast<int>(table.size())) return table.back();
    const double f = u - k;
    return table[k] * (1.0 - f) + table[k + 1] * f;
}

double Kilonova::luminosity(double t) const {
    return sample(m_L, t);
}
double Kilonova::temperature(double t) const {
    return sample(m_T, t);
}
double Kilonova::photosphere(double t) const {
    return sample(m_R, t);
}

// --- The jet --------------------------------------------------------------------------------------

double breakoutTime(const Jet& jet, double vmax) {
    if (jet.headSpeed <= vmax) return 1e30;
    return jet.headSpeed * jet.launch / (jet.headSpeed - vmax);
}

double dopplerFactor(double lorentz, double angle) {
    const double beta = std::sqrt(std::max(1.0 - 1.0 / (lorentz * lorentz), 0.0));
    return 1.0 / (lorentz * (1.0 - beta * std::cos(angle)));
}

// --- Time --------------------------------------------------------------------------------------

double tempo(const Tempo& p, double t, double collapse, double burst) {
    double rate;
    if (t < 0.0) {
        rate = p.inspiralRate * -t;
    } else {
        // Smoothly from the slow pace of the hypermassive star to the faster one after the
        // collapse, over a factor of four in time.
        const double s = std::clamp(std::log(std::max(t, 1e-12) / collapse) / std::log(4.0), 0.0, 1.0);
        const double k = p.afterRate + (p.lateRate - p.afterRate) * s * s * (3.0 - 2.0 * s);
        rate = k * t;
        // Hold for the gamma-ray burst: slow down across a few tenths of a decade around it.
        if (burst > 0.0 && t > 0.0) {
            const double d = std::log(t / burst) / 0.18;
            rate /= 1.0 + p.burstHold * 4.0 * std::exp(-d * d);
        }
    }
    return std::max(rate, p.mergerRate);
}

}  // namespace merger
