#pragma once

// What stars and planets look like, from the physics that shapes their surfaces: limb
// darkening, the size and lifetime of convection cells, starspots, the corona and the
// sky of an Earth-like planet. The ray tracer draws them (shaders/surfaces.glsl and
// shaders/central_surface.glsl); this is the tested CPU side that sets their scales.

#include "central.h"

#include <glm/glm.hpp>

#include <vector>

namespace surface {

// --- Limb darkening -------------------------------------------------------------------

// Eddington-Barbier with a grey atmosphere: light leaving at mu = cos(angle from the
// normal) comes from optical depth ~mu, where T^4 = (3/4) Teff^4 (tau + 2/3). Returns
// T(mu) / Teff: the limb shows cooler, higher layers. (Integrated over the spectrum this
// is Eddington's I(mu) / I(1) = (2 + 3 mu) / 5.)
double limbTemperatureFactor(double mu);
// Planck's spectral radiance B_lambda(T), W / (m^2 sr m).
double planck(double lambdaNm, double temperature);
// I_lambda(mu) / I_lambda(1) for a grey atmosphere: darker and redder toward the limb.
double limbDarkening(double teff, double mu, double lambdaNm);

// --- Convection ----------------------------------------------------------------------

// Newtonian surface gravity of a body of mass M (solar masses) and radius R (km), cm/s^2.
double surfaceGravityCgs(double massSolar, double radiusKm);
// Pressure scale height of a photosphere, k T / (mu m_H g), with mean molecular weight
// 1.3 (neutral hydrogen and helium), in km.
double pressureScaleHeightKm(double temperature, double gravityCgs);
// Adiabatic sound speed sqrt(gamma k T / (mu m_H)), gamma = 5/3, in km/s.
double soundSpeedKmS(double temperature);
// Granules are ~10 pressure scale heights across (~1,300 km on the Sun) and last as long
// as it takes gas rising at ~1/5 of the sound speed to cross one (~14 minutes on the
// Sun); on a red giant a handful of cells cover the star and live for months.
double granuleSizeKm(double temperature, double gravityCgs);
double granuleLifetimeSeconds(double temperature, double gravityCgs);

// --- Starspots ---------------------------------------------------------------------

// Umbral temperature drop below the photosphere, in K: ~1,800 K on the Sun, only a few
// hundred on cool red dwarfs (the trend of Berdyugina 2005, Living Rev. Solar Phys. 2, 8).
double spotTemperatureDrop(double temperature);

struct Spot {
    glm::vec3 direction;  // in the star's rotating frame (z along the spin axis)
    float radius;         // angular radius, radians
    float umbra;          // umbra radius / spot radius (1 for a pore: no penumbra)
};
struct ActiveRegion {
    glm::vec3 direction;
    float radius;  // angular radius of its plage, radians
};
constexpr int kMaxSpots = 32, kMaxRegions = 8;

enum class SpotStyle {
    Solar,      // bipolar groups in two belts below ~35 degrees (the butterfly diagram)
    RedDwarf,   // fewer, larger spots at all latitudes, often one near a pole
    Giant,      // a few huge spots
};
struct SpotPattern {
    std::vector<Spot> spots;
    std::vector<ActiveRegion> regions;
};
// activity: 0 (quiet) to 1 (solar maximum or an active star). Spots of a group: the
// leader, closer to the equator and ahead in the direction of rotation (Joy's law: the
// group is tilted by ~half its latitude), the follower behind it, and small pores.
SpotPattern makeSpots(SpotStyle style, double activity, unsigned seed);

// --- Corona ------------------------------------------------------------------------

// Electron density of the solar corona (Baumbach and Allen), in 1e8 per cm^3, at
// rho = r / R: 2.99 rho^-16 + 1.55 rho^-6 + 0.036 rho^-1.5.
double coronaElectronDensity(double rho);
// The dilution of the photosphere's light at rho: (1 - sqrt(1 - 1/rho^2)) / 2.
double dilution(double rho);
// Thomson-scattered light along a line of sight passing rho from the center: the
// integral of n_e W along it (in 1e8 cm^-3 R).
double coronaColumn(double rho);
// Baumbach's measured K-corona brightness, relative to the mean brightness of the disk:
// 1e-6 (0.0532 rho^-2.5 + 1.425 rho^-7 + 2.565 rho^-17).
double baumbachBrightness(double rho);
// Converts the column into brightness relative to the disk center (Eddington limb
// darkening: the mean disk is 0.8 of its center), for the real Sun.
double coronaScaleReal();

// --- The appearance of the central object --------------------------------------------

enum class Look { None, Photosphere, BrownDwarf, WhiteDwarf, NeutronStar, Earthlike, GasGiant };

struct Appearance {
    Look look = Look::None;
    double granuleFrequency = 0.0;  // convection cells per radian of the surface
    double granuleContrast = 0.0;   // rms temperature contrast of the cells
    double granuleLane = 0.15;      // width of the dark lanes, relative to a cell
    double granuleLifetime = 0.0;   // seconds
    double spotDrop = 0.0;          // umbral temperature drop / surface temperature
    SpotStyle spotStyle = SpotStyle::Solar;
    double activity = 0.0;
    bool corona = false;
    bool atmosphere = false;        // an Earth-like sky
    double flareInterval = 0.0;     // mean seconds between white-light flares (0: none)
};
Appearance appearance(const central::Object& object);

}  // namespace surface
