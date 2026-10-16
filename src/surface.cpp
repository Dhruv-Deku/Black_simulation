#include "surface.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace surface {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kBoltzmann = 1.380649e-16;     // erg / K
constexpr double kHydrogenMass = 1.6735575e-24; // g
constexpr double kGravitationalConstant = 6.6743e-8;  // cgs
constexpr double kSolarMassGrams = 1.98847e33;
constexpr double kMeanMolecularWeight = 1.3;    // neutral hydrogen and helium

glm::vec3 direction(double latitude, double longitude) {
    return glm::vec3(static_cast<float>(std::cos(latitude) * std::cos(longitude)),
                     static_cast<float>(std::cos(latitude) * std::sin(longitude)),
                     static_cast<float>(std::sin(latitude)));
}

}  // namespace

// --- Limb darkening -------------------------------------------------------------------

double limbTemperatureFactor(double mu) {
    return std::pow(0.75 * (std::clamp(mu, 0.0, 1.0) + 2.0 / 3.0), 0.25);
}

double planck(double lambdaNm, double temperature) {
    constexpr double h = 6.62607015e-34, c = 2.99792458e8, k = 1.380649e-23;
    const double lambda = lambdaNm * 1e-9;
    return 2.0 * h * c * c / std::pow(lambda, 5.0) / std::expm1(h * c / (lambda * k * temperature));
}

double limbDarkening(double teff, double mu, double lambdaNm) {
    return planck(lambdaNm, teff * limbTemperatureFactor(mu)) / planck(lambdaNm, teff * limbTemperatureFactor(1.0));
}

// --- Convection ----------------------------------------------------------------------

double surfaceGravityCgs(double massSolar, double radiusKm) {
    const double r = radiusKm * 1e5;
    return kGravitationalConstant * massSolar * kSolarMassGrams / (r * r);
}

double pressureScaleHeightKm(double temperature, double gravityCgs) {
    return kBoltzmann * temperature / (kMeanMolecularWeight * kHydrogenMass * gravityCgs) / 1e5;
}

double soundSpeedKmS(double temperature) {
    return std::sqrt(5.0 / 3.0 * kBoltzmann * temperature / (kMeanMolecularWeight * kHydrogenMass)) / 1e5;
}

double granuleSizeKm(double temperature, double gravityCgs) {
    return 10.0 * pressureScaleHeightKm(temperature, gravityCgs);
}

double granuleLifetimeSeconds(double temperature, double gravityCgs) {
    return granuleSizeKm(temperature, gravityCgs) / (0.2 * soundSpeedKmS(temperature));
}

// --- Starspots ---------------------------------------------------------------------

double spotTemperatureDrop(double temperature) {
    // Through Berdyugina's compilation: ~2,000 K for G stars, a few hundred for M dwarfs.
    return 1850.0 * std::exp((temperature - 5772.0) / 1200.0);
}

SpotPattern makeSpots(SpotStyle style, double activity, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);
    const double deg = kPi / 180.0;
    SpotPattern out;
    auto addSpot = [&](double lat, double lon, double radius, double umbra) {
        if (static_cast<int>(out.spots.size()) >= kMaxSpots) return;
        out.spots.push_back({direction(lat, lon), static_cast<float>(radius), static_cast<float>(umbra)});
    };
    auto addRegion = [&](double lat, double lon, double radius) {
        if (static_cast<int>(out.regions.size()) >= kMaxRegions) return;
        out.regions.push_back({direction(lat, lon), static_cast<float>(radius)});
    };
    activity = std::clamp(activity, 0.0, 1.0);

    if (style == SpotStyle::Solar) {
        // Groups emerge in two belts that drift from ~30 to ~8 degrees over the cycle
        // (the butterfly diagram); mid-cycle they sit near 15 degrees.
        const int groups = static_cast<int>(std::lround(1.0 + 8.0 * activity));
        for (int g = 0; g < groups; ++g) {
            const double hemisphere = uni(rng) < 0.5 ? -1.0 : 1.0;
            const double lat = hemisphere * std::clamp(15.0 + 6.0 * normal(rng), 5.0, 38.0) * deg;
            // The first, big group faces the default camera (longitude 0) at the start.
            const double lon = g == 0 ? 0.35 : 2.0 * kPi * uni(rng);
            // Big spots are ~10,000-30,000 km across (the Sun's radius is 696,000 km).
            // The first group of the set is a big one, like those seen near solar maximum.
            const double leader = g == 0 ? 0.03 + 0.012 * uni(rng)
                                         : std::clamp(std::exp(std::log(0.014) + 0.45 * normal(rng)), 0.005, 0.045);
            const double separation = 5.0 * leader + (0.03 + 0.06 * uni(rng));
            // Joy's law: the leader sits nearer the equator, the group tilted by ~lat/2.
            const double tilt = 0.5 * std::abs(lat);
            const double leadLat = lat - hemisphere * 0.5 * separation * std::sin(tilt);
            const double leadLon = lon + 0.5 * separation * std::cos(tilt) / std::cos(lat);
            const double followLat = lat + hemisphere * 0.5 * separation * std::sin(tilt);
            const double followLon = lon - 0.5 * separation * std::cos(tilt) / std::cos(lat);
            addSpot(leadLat, leadLon, leader, 0.38 + 0.08 * uni(rng));
            if (uni(rng) < 0.8) addSpot(followLat, followLon, leader * (0.45 + 0.35 * uni(rng)), 0.4);
            // Pores: small spots without a penumbra, scattered between them.
            const int pores = 1 + static_cast<int>(3.0 * uni(rng));
            for (int k = 0; k < pores; ++k) {
                const double f = uni(rng);
                addSpot(std::lerp(leadLat, followLat, f) + 0.25 * separation * normal(rng),
                        std::lerp(leadLon, followLon, f) + 0.25 * separation * normal(rng), 0.0025 + 0.003 * uni(rng), 1.0);
            }
            addRegion(lat, lon, 0.6 * separation + 2.5 * leader);
        }
    } else if (style == SpotStyle::RedDwarf) {
        // Fully convective red dwarfs: fewer but far larger spots, at any latitude, and
        // often one over a pole (seen in Doppler images).
        const int count = static_cast<int>(std::lround(2.0 + 4.0 * activity));
        for (int k = 0; k < count; ++k) {
            const double lat = std::asin(std::clamp(1.9 * uni(rng) - 0.95, -0.95, 0.95));
            const double lon = k == 0 ? 0.4 : 2.0 * kPi * uni(rng);  // one in view at the start
            const double radius = 0.05 + 0.1 * uni(rng);
            addSpot(lat, lon, radius, 0.55 + 0.2 * uni(rng));
            addRegion(lat, lon, 1.8 * radius);
        }
        if (uni(rng) < 0.6 + 0.3 * activity) {
            const double lat = (uni(rng) < 0.5 ? -1.0 : 1.0) * (75.0 + 13.0 * uni(rng)) * deg;
            addSpot(lat, 2.0 * kPi * uni(rng), 0.14 + 0.1 * uni(rng), 0.6);
        }
    } else {
        // Giants: a few huge spots.
        const int count = static_cast<int>(std::lround(2.0 + 2.0 * activity));
        for (int k = 0; k < count; ++k) {
            const double lat = std::clamp(35.0 * normal(rng), -65.0, 65.0) * deg;
            const double lon = k == 0 ? 0.5 : 2.0 * kPi * uni(rng);
            const double radius = 0.12 + 0.15 * uni(rng);
            addSpot(lat, lon, radius, 0.5);
            addRegion(lat, lon, 1.6 * radius);
        }
    }
    return out;
}

// --- Corona ------------------------------------------------------------------------

double coronaElectronDensity(double rho) {
    return 2.99 * std::pow(rho, -16.0) + 1.55 * std::pow(rho, -6.0) + 0.036 * std::pow(rho, -1.5);
}

double dilution(double rho) {
    return 0.5 * (1.0 - std::sqrt(std::max(1.0 - 1.0 / (rho * rho), 0.0)));
}

double coronaColumn(double rho) {
    // 2 * integral over z >= 0 of n_e W at r = sqrt(rho^2 + z^2); with z = rho tan(psi)
    // the integrand is smooth on [0, pi/2) and falls off as cos(psi)^1.5.
    const int n = 4000;
    const double h = 0.5 * kPi / n;
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double psi = std::min(i * h, 0.5 * kPi - 1e-9);
        const double c = std::cos(psi);
        const double r = rho / c;
        const double f = coronaElectronDensity(r) * dilution(r) * rho / (c * c);
        const double w = (i == 0 || i == n) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum += w * f;
    }
    return 2.0 * sum * h / 3.0;
}

double baumbachBrightness(double rho) {
    return 1e-6 * (0.0532 * std::pow(rho, -2.5) + 1.425 * std::pow(rho, -7.0) + 2.565 * std::pow(rho, -17.0));
}

double coronaScaleReal() {
    // Fix the scale at 1.1 R, where both the density and the brightness are best measured.
    return 0.8 * baumbachBrightness(1.1) / coronaColumn(1.1);
}

// --- Appearance -----------------------------------------------------------------------

Appearance appearance(const central::Object& o) {
    using central::Kind;
    Appearance a;
    switch (o.kind) {
        case Kind::BlackHole: a.look = Look::None; return a;
        case Kind::NeutronStar: a.look = Look::NeutronStar; return a;
        case Kind::WhiteDwarf: a.look = Look::WhiteDwarf; return a;
        case Kind::RockyPlanet:
            a.look = Look::Earthlike;
            a.atmosphere = true;
            return a;
        case Kind::GasGiant: a.look = Look::GasGiant; return a;
        case Kind::Star: break;
    }
    const double T = o.surfaceTemperature;
    if (T < 2500.0) {
        // Brown dwarfs (and the coolest stars): clouds instead of granulation and spots.
        a.look = Look::BrownDwarf;
        return a;
    }
    a.look = Look::Photosphere;
    const double radiusKm = central::lengthKm(o.radius, o.massSolar);
    const double g = surfaceGravityCgs(o.massSolar, radiusKm);
    a.granuleFrequency = radiusKm / granuleSizeKm(T, g);
    a.granuleLifetime = granuleLifetimeSeconds(T, g);
    const double logg = std::log10(g);
    const bool giant = logg < 3.5;
    // rms temperature contrast: ~3% on the Sun (13% in intensity at 550 nm); giants'
    // few huge cells are more vigorous.
    a.granuleContrast = giant ? 0.05 : 0.029;
    a.granuleLane = giant ? 0.6 : 0.15;
    a.spotDrop = std::min(spotTemperatureDrop(T) / T, 0.6);
    if (giant) {
        a.spotStyle = SpotStyle::Giant;
        a.activity = 0.4;
    } else if (T < 3900.0) {
        a.spotStyle = SpotStyle::RedDwarf;
        a.activity = 0.8;
        a.flareInterval = 2.0 * 3600.0;  // an active red dwarf like Proxima flares every few hours
    } else {
        a.spotStyle = SpotStyle::Solar;
        a.activity = 0.85;  // near the maximum of solar cycle 25 (2024-26)
    }
    a.corona = !giant;  // dwarfs have hot coronae; red giants' are weak
    return a;
}

}  // namespace surface
