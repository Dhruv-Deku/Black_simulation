#include "tde.h"

#include "central.h"

#include <algorithm>
#include <cmath>

namespace tde {

namespace {

constexpr double kTwoPi = 2.0 * kerr::kPi;

// cgs constants for the physical scales.
constexpr double kG = 6.6743e-8;             // cm^3 g^-1 s^-2
constexpr double kSolarMassGrams = 1.98847e33;
constexpr double kYearSeconds = 3.15576e7;   // Julian year

// R(r) = (dr/dtau)^2 for mu = 1, M = 1, written as the quartic
// (E^2 - 1) r^4 + 2 r^3 + (a^2 (E^2 - 1) - L^2 - Q) r^2 + 2 ((L - aE)^2 + Q) r - a^2 Q.
double radial(double a, double E, double L, double Q, double r) {
    const double x = L - a * E;
    return (((E * E - 1.0) * r + 2.0) * r + (a * a * (E * E - 1.0) - L * L - Q)) * r * r +
           2.0 * (x * x + Q) * r - a * a * Q;
}

// Bisection for the zero of R between lo (R >= 0) and hi (R < 0).
double bisect(double a, double E, double L, double Q, double lo, double hi) {
    for (int k = 0; k < 90; ++k) {
        const double mid = 0.5 * (lo + hi);
        (radial(a, E, L, Q, mid) >= 0.0 ? lo : hi) = mid;
        if (std::abs(hi - lo) < 1e-12 * hi) break;
    }
    return 0.5 * (lo + hi);
}

}  // namespace

Scales classicalScales(double M, double m, double rStar) {
    Scales s;
    s.tidalRadius = rStar * std::cbrt(M / m);
    s.energySpread = M * rStar / (s.tidalRadius * s.tidalRadius);
    s.aMin = M / (2.0 * s.energySpread);
    s.tMin = kTwoPi * std::sqrt(s.aMin * s.aMin * s.aMin / M);
    s.peakRate = m / (3.0 * s.tMin);
    return s;
}

Scales physicalScales(double holeSolar, double mStar, double rStarKm) {
    // In cgs with G: lengths in cm, then reported in km; times in s; rates in Suns a year.
    const double GM = kG * holeSolar * kSolarMassGrams;
    const double r = rStarKm * 1e5;
    Scales s;
    const double rt = r * std::cbrt(holeSolar / mStar);
    s.energySpread = GM * r / (rt * rt);
    const double aMin = GM / (2.0 * s.energySpread);
    s.tMin = kTwoPi * std::sqrt(aMin * aMin * aMin / GM);
    s.peakRate = mStar / (3.0 * s.tMin) * kYearSeconds;
    s.tidalRadius = rt * 1e-5;
    s.aMin = aMin * 1e-5;
    return s;
}

double keplerPeriod(double M, double e) {
    return kTwoPi * M / std::pow(-2.0 * e, 1.5);
}

RadialOrbit radialOrbit(double a, const kerr::Constants& c, const kerr::State& s) {
    RadialOrbit out;
    const double E = c.E, L = c.L, Q = c.Q;
    const double rh = kerr::horizonRadius(a);
    double r0 = s.r;
    // Parcels handed over in single precision can sit a hair outside the allowed region
    // (R slightly negative): step back into it.
    if (radial(a, E, L, Q, r0) < 0.0) {
        for (int k = 1; k <= 40; ++k) {
            const double d = r0 * 5e-4 * k;
            if (radial(a, E, L, Q, r0 + d) >= 0.0) {
                r0 += d;
                break;
            }
            if (r0 - d > rh && radial(a, E, L, Q, r0 - d) >= 0.0) {
                r0 -= d;
                break;
            }
        }
    }

    // The nearest turning point inward: scan down to the horizon, then bisect. R(r+) =
    // P(r+)^2 >= 0, so an inner zero shows up as a sign change on the way down.
    bool hasInner = false;
    {
        const int n = 160;
        const double ratio = std::pow((rh * (1.0 + 1e-6)) / r0, 1.0 / n);
        double prev = r0;
        for (int i = 1; i <= n; ++i) {
            const double r = r0 * std::pow(ratio, i);
            if (radial(a, E, L, Q, r) < 0.0) {
                out.rp = bisect(a, E, L, Q, prev, r);
                hasInner = true;
                break;
            }
            prev = r;
        }
    }
    const bool inbound = s.vr < 0.0;
    if (E >= 1.0) {
        out.kind = (!hasInner && inbound) ? RadialOrbit::Plunges : RadialOrbit::Unbound;
        return out;
    }
    if (!hasInner) {
        // Bound, but nothing stops it short of the horizon: it falls in on this pass
        // (or the next, if it is heading out now).
        out.kind = RadialOrbit::Plunges;
        return out;
    }

    // The outer turning point: R < 0 beyond 2 / (1 - E^2), where the r^3 and r^4 terms
    // already make it negative; only one zero lies above a radius inside the orbit.
    const double e2 = 1.0 - E * E;
    double hi = std::max(2.0 / e2 + 4.0, r0 * 1.01);
    for (int k = 0; k < 60 && radial(a, E, L, Q, hi) >= 0.0; ++k) hi *= 2.0;
    out.ra = bisect(a, E, L, Q, r0, hi);
    out.kind = RadialOrbit::Bound;

    // With rp and ra known, R = (1 - E^2)(ra - r)(r - rp)(r^2 - s r + q): the sum of the
    // four roots is 2 / (1 - E^2) and their product a^2 Q / (1 - E^2).
    const double rp = out.rp, ra = out.ra;
    const double sum = 2.0 / e2 - ra - rp;
    const double prod = a * a * Q / (e2 * ra * rp);
    // Average of sin^2(theta) over the polar oscillation, whose turning points are the
    // root z of a^2 (1 - E^2) z^2 - (Q + L^2 + a^2 (1 - E^2)) z + Q = 0, z = cos^2(theta).
    double zMax = 0.0;
    {
        const double A = a * a * e2, B = Q + L * L + A;
        zMax = A > 1e-12 ? (B - std::sqrt(std::max(B * B - 4.0 * A * Q, 0.0))) / (2.0 * A) : Q / std::max(B, 1e-300);
    }
    const double sin2 = 1.0 - 0.5 * std::clamp(zMax, 0.0, 1.0);
    const double rm = 0.5 * (ra + rp), h = 0.5 * (ra - rp);
    // dt/dchi: (dt/dtau) / sqrt((1 - E^2)(r^2 - s r + q)), with Boyer-Lindquist
    // dt/dtau = (r^2 + a^2) P / Delta + a L - a^2 E <sin^2 theta> (the same between
    // pericenters in Kerr-Schild time, which differs by a function of r alone).
    auto rate = [&](double chi) {
        const double r = rm - h * std::cos(chi);
        const double rra = r * r + a * a;
        const double P = E * rra - a * L;
        const double delta = r * r - 2.0 * r + a * a;
        const double dtdtau = rra * P / delta + a * L - a * a * E * sin2;
        const double quad = std::max(r * r - sum * r + prod, 1e-12 * r * r);
        return dtdtau / std::sqrt(e2 * quad);
    };
    // Simpson's rule in chi from 0 to x.
    auto integral = [&](double x, int n) {
        if (x <= 0.0) return 0.0;
        const double step = x / n;
        double total = rate(0.0) + rate(x);
        for (int i = 1; i < n; ++i) total += (i % 2 ? 4.0 : 2.0) * rate(i * step);
        return total * step / 3.0;
    };
    out.period = 2.0 * integral(kerr::kPi, 96);
    const double chi0 = std::acos(std::clamp((rm - r0) / h, -1.0, 1.0));
    const double sincePericenter = integral(chi0, 64);
    // In Boyer-Lindquist time; the simulation's Kerr-Schild clock runs ahead by
    // timeShift(r), which differs between r0 and the pericenter.
    out.toPericenter = (inbound ? sincePericenter : out.period - sincePericenter) + kerr::timeShift(a, rp) -
                       kerr::timeShift(a, r0);
    return out;
}

void Fallback::build(const std::vector<double>& returnTimes, const std::vector<double>& weights, double totalWeight) {
    m_rate.clear();
    m_cumulative.clear();
    boundFraction = firstReturn = peakTime = peakRate = lateSlope = 0.0;
    if (returnTimes.empty() || totalWeight <= 0.0) return;
    double sum = 0.0;
    double tMin = 1e300, tMax = 0.0;
    for (size_t i = 0; i < returnTimes.size(); ++i) {
        sum += weights[i];
        tMin = std::min(tMin, returnTimes[i]);
        tMax = std::max(tMax, returnTimes[i]);
    }
    boundFraction = sum / totalWeight;
    firstReturn = std::max(tMin, 1e-6);

    // A grid in ln t from a quarter of the first return out to the last (at most 3000
    // first returns: the tail beyond falls off as t^(-5/3) and is extrapolated).
    const int n = 900;
    const double t0 = 0.25 * firstReturn;
    const double t1 = std::clamp(tMax * 1.5, 4.0 * firstReturn, 3000.0 * firstReturn);
    m_logT0 = std::log(t0);
    m_dLogT = (std::log(t1) - m_logT0) / (n - 1);
    std::vector<double> mass(n, 0.0);
    for (size_t i = 0; i < returnTimes.size(); ++i) {
        const double x = (std::log(std::max(returnTimes[i], t0)) - m_logT0) / m_dLogT;
        if (x > n - 1) continue;  // beyond the table: part of the extrapolated tail
        const int j = std::clamp(static_cast<int>(std::lround(x)), 0, n - 1);
        mass[j] += weights[i] / totalWeight;
    }
    // Smooth over 4% in time: a finite number of parcels gives a noisy histogram.
    const double sigma = 0.04 / m_dLogT;
    const int reach = static_cast<int>(std::ceil(4.0 * sigma));
    std::vector<double> kernel(2 * reach + 1);
    double norm = 0.0;
    for (int k = -reach; k <= reach; ++k) norm += kernel[k + reach] = std::exp(-0.5 * k * k / (sigma * sigma));
    for (double& w : kernel) w /= norm;
    std::vector<double> smooth(n, 0.0);
    for (int j = 0; j < n; ++j) {
        if (mass[j] == 0.0) continue;
        for (int k = -reach; k <= reach; ++k) {
            const int i = j + k;
            if (i >= 0 && i < n) smooth[i] += mass[j] * kernel[k + reach];
        }
    }
    m_rate.resize(n);
    m_cumulative.resize(n);
    double cumulative = 0.0;
    for (int j = 0; j < n; ++j) {
        const double t = std::exp(m_logT0 + j * m_dLogT);
        m_rate[j] = smooth[j] / (m_dLogT * t);  // dN/dt = (dN/d ln t) / t
        cumulative += smooth[j];
        m_cumulative[j] = cumulative;
        if (m_rate[j] > peakRate) {
            peakRate = m_rate[j];
            peakTime = t;
        }
    }
    // The decline: a least-squares slope of ln(rate) against ln(t) from 3 to 20 peak times.
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    int count = 0;
    for (int j = 0; j < n; ++j) {
        const double t = std::exp(m_logT0 + j * m_dLogT);
        if (t < 3.0 * peakTime || t > 20.0 * peakTime || m_rate[j] <= 0.0) continue;
        const double x = std::log(t), y = std::log(m_rate[j]);
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        ++count;
    }
    if (count > 5) lateSlope = (count * sxy - sx * sy) / (count * sxx - sx * sx);
}

double Fallback::rate(double t) const {
    if (m_rate.empty() || t <= 0.0) return 0.0;
    const double x = (std::log(t) - m_logT0) / m_dLogT;
    if (x <= 0.0) return x > -1.0 ? m_rate.front() * (1.0 + x) : 0.0;
    const int last = static_cast<int>(m_rate.size()) - 1;
    if (x >= last) {
        const double tLast = std::exp(m_logT0 + last * m_dLogT);
        return m_rate[last] * std::pow(t / tLast, -5.0 / 3.0);
    }
    const int j = static_cast<int>(x);
    const double f = x - j;
    return m_rate[j] * (1.0 - f) + m_rate[j + 1] * f;
}

double Fallback::returned(double t) const {
    if (m_cumulative.empty() || t <= 0.0) return 0.0;
    const double x = (std::log(t) - m_logT0) / m_dLogT;
    if (x <= 0.0) return 0.0;
    const int last = static_cast<int>(m_cumulative.size()) - 1;
    if (x >= last) {
        // The t^(-5/3) tail: integral of r_L (t / t_L)^(-5/3) dt from t_L to t.
        const double tLast = std::exp(m_logT0 + last * m_dLogT);
        return m_cumulative[last] + 1.5 * m_rate[last] * tLast * (1.0 - std::pow(t / tLast, -2.0 / 3.0));
    }
    const int j = static_cast<int>(x);
    const double f = x - j;
    return m_cumulative[j] * (1.0 - f) + m_cumulative[j + 1] * f;
}

double eddingtonLuminosity(double massSolar) {
    // 4 pi G M m_p c / sigma_T.
    return 1.2572e38 * massSolar;
}

double escapingLuminosity(double lacc) {
    return lacc <= 1.0 ? std::max(lacc, 0.0) : 1.0 + std::log(lacc);
}

double windSpeed(double M, double vInf, double R) {
    return std::sqrt(vInf * vInf + 2.0 * M / R);
}

double windTravelTime(double M, double vInf, double R0, double R) {
    const double k = vInf * vInf;
    auto T = [&](double x) {
        return std::sqrt(x * (k * x + 2.0 * M)) / k -
               2.0 * M / (k * std::sqrt(k)) * std::log(std::sqrt(k * x) + std::sqrt(k * x + 2.0 * M));
    };
    return T(R) - T(R0);
}

double windFront(double M, double vInf, double R0, double t) {
    if (t <= 0.0) return R0;
    // The front moves no faster than it started and no slower than vInf.
    double lo = R0, hi = R0 + windSpeed(M, vInf, R0) * t + 1.0;
    for (int k = 0; k < 80; ++k) {
        const double mid = 0.5 * (lo + hi);
        (windTravelTime(M, vInf, R0, mid) < t ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double openingAngle(double late, double first, double formed) {
    return late + (first - late) * std::exp(-std::max(formed, 0.0) / 0.2);
}

Event analyze(double a, double time, const std::vector<ParticleInit>& debris, int maxSamples) {
    Event ev;
    ev.time = time;
    if (debris.empty()) return ev;
    const size_t stride = std::max<size_t>(1, debris.size() / std::max(maxSamples, 1));
    std::vector<double> times, weights;
    int total = 0, unbound = 0, plunging = 0;
    for (size_t i = 0; i < debris.size(); i += stride) {
        const ParticleInit& p = debris[i];
        const kerr::Constants c{p.E, p.L, p.Q, 1.0};
        const kerr::State s{time, p.r, p.theta, p.phi, p.vr, p.vtheta};
        const RadialOrbit orbit = radialOrbit(a, c, s);
        ++total;
        switch (orbit.kind) {
            case RadialOrbit::Unbound: ++unbound; break;
            case RadialOrbit::Plunges:
                // Swallowed within about a free-fall time: accreted at once.
                ++plunging;
                times.push_back(0.5 * std::pow(p.r, 1.5));
                weights.push_back(1.0);
                break;
            case RadialOrbit::Bound: {
                // Debris is torn off near pericenter: the passage it is making does not
                // count as a return; the next one does.
                const bool passing = p.vr < 0.0 && p.r < 3.0 * orbit.rp;
                times.push_back(passing ? orbit.toPericenter + orbit.period : orbit.toPericenter);
                weights.push_back(1.0);
                break;
            }
        }
    }
    ev.unboundFraction = static_cast<double>(unbound) / total;
    ev.plungingFraction = static_cast<double>(plunging) / total;
    ev.fallback.build(times, weights, total);
    return ev;
}

Aftermath::Aftermath(const std::vector<Event>& events, const OutflowParams& params, double holeSolar, double spin)
    : m_events(events), m_params(params) {
    double strongest = 0.0;
    for (size_t i = 0; i < m_events.size(); ++i) {
        const Event& e = m_events[i];
        if (e.fallback.empty()) continue;
        const double peak = e.fallback.peakRate * e.starMass;
        if (peak > strongest) {
            strongest = peak;
            m_reference = static_cast<int>(i);
        }
    }
    if (m_reference < 0) return;
    // The summed rate peaks near one of the events' own peaks.
    auto sum = [&](double t) {
        double total = 0.0;
        for (const Event& e : m_events) total += e.fallback.rate(t - e.time) * e.starMass;
        return total;
    };
    double first = 1e300;
    for (const Event& e : m_events) {
        if (e.fallback.empty()) continue;
        first = std::min(first, e.time + 0.25 * e.fallback.firstReturn);
        for (int k = 0; k <= 200; ++k) {
            const double t = e.time + e.fallback.peakTime * (0.6 + 0.004 * k);
            const double r = sum(t);
            if (r > m_peakMass) {
                m_peakMass = r;
                m_peakTime = t;
            }
        }
    }
    // When the returning debris starts to arrive: the rate first passes 1% of its peak.
    m_start = m_peakTime;
    for (int k = 0; k <= 4000; ++k) {
        const double t = first + (m_peakTime - first) * k / 4000.0;
        if (sum(t) > 0.01 * m_peakMass) {
            m_start = t;
            break;
        }
    }
    // The flare's scale: a Sun-like star torn apart by this hole returns at most
    // physicalScales().peakRate; accreting the share that is not blown out releases the
    // hole's efficiency (1 - E_isco) times its rest energy.
    const Scales sun = physicalScales(holeSolar, 1.0, central::kSunRadiusKm);
    const double mdot = sun.peakRate * kSolarMassGrams / kYearSeconds;  // g/s
    constexpr double kC2 = 8.987551787e20;                                // c^2, erg/g
    m_peakAccretion = central::blackHoleEfficiency(spin) * (1.0 - m_params.fraction) * mdot * kC2 /
                      eddingtonLuminosity(holeSolar);
}

double Aftermath::level(double t) const {
    if (m_peakMass <= 0.0) return 0.0;
    double total = 0.0;
    for (const Event& e : m_events) total += e.fallback.rate(t - e.time) * e.starMass;
    return total / m_peakMass;
}

double Aftermath::returned(double t) const {
    if (!active()) return 0.0;
    double total = 0.0;
    for (const Event& e : m_events) total += e.fallback.returned(t - e.time) * e.starMass;
    return total / reference().starMass;
}

double Aftermath::opening(double t) const {
    const double first = std::min(m_params.openingLate * m_params.widening, 1.35);
    const double bound = std::max(reference().fallback.boundFraction, 1e-6);
    return openingAngle(m_params.openingLate, std::max(first, m_params.openingLate), returned(t) / bound);
}

double Aftermath::accretionRatio(double t) const {
    return m_peakAccretion * level(t);
}

double Aftermath::eddingtonRatio(double t) const {
    return escapingLuminosity(accretionRatio(t));
}

double Aftermath::luminosity(double t) const {
    const double peak = escapingLuminosity(m_peakAccretion);
    return peak > 0.0 ? eddingtonRatio(t) / peak : 0.0;
}

double Aftermath::launchRadius() const {
    return active() ? std::max(reference().circularization, 10.0) : 10.0;
}

double Aftermath::front(double t) const {
    return windFront(1.0, m_params.speed, launchRadius(), t - m_start);
}

}  // namespace tde
