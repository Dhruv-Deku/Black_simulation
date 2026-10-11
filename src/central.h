#pragma once

// The object at the center: a black hole, a neutron star, a white dwarf, a star or a
// planet.
//
// Units stay G = c = M = 1 with M the central mass, so lengths are in GM/c^2 and times
// in GM/c^3. The physical mass (in solar masses) only sets what those units are in
// kilometres and seconds.
//
// Outside any spherical, non-rotating body spacetime is exactly Schwarzschild
// (Birkhoff's theorem): for the same mass, orbits and light bending outside the
// surface are identical for a black hole, a neutron star and the Sun. Only how close
// one can get changes. A rotating body is modelled by the Kerr exterior with the
// body's mass and angular momentum, cut off at its surface: right in M and J, but
// with the Kerr quadrupole (a neutron star's real quadrupole is 2-10 times larger, a
// small effect for spins below ~0.4).

#include <string>
#include <vector>

namespace central {

// --- Physical scale --------------------------------------------------------------

constexpr double kKmPerSolarMass = 1.476625;        // G M_sun / c^2
constexpr double kSecondsPerSolarMass = 4.925491e-6;  // G M_sun / c^3
constexpr double kSunRadiusKm = 695700.0;
constexpr double kAuKm = 1.495978707e8;
constexpr double kLightYearKm = 9.4607304725808e12;
constexpr double kSpeedOfLightKmS = 299792.458;
constexpr double kChandrasekharMass = 1.44;           // M_sun (mu_e = 2)
constexpr double kBuchdahlRadius = 9.0 / 4.0;         // no static star is smaller (in M)

double lengthKm(double lengthM, double massSolar);
double timeSeconds(double timeM, double massSolar);
double accelerationSI(double accelerationGeometric, double massSolar);  // m/s^2
// "1,234 km", "3.2 AU", "4.1 ly"; "12 us", "3.4 s", "2.1 h", "5.2 yr".
std::string formatLength(double lengthM, double massSolar);
std::string formatTime(double timeM, double massSolar);

// --- The object --------------------------------------------------------------------

enum class Kind {
    BlackHole,
    NeutronStar,
    WhiteDwarf,
    Star,         // main sequence, red dwarf, brown dwarf, red giant
    RockyPlanet,
    GasGiant,
};
const char* kindName(Kind kind);

struct Object {
    Kind kind = Kind::BlackHole;
    std::string name = "Black hole";
    double massSolar = 4.3e6;  // physical mass in M_sun (sets the units only)
    double radius = 0.0;       // surface radius R in M (unused for a black hole)
    double chi = 0.94;         // spin c J / (G M^2), physical (planets exceed 1)
    double inertia = 0.35;     // moment of inertia I / (M R^2)
    double surfaceTemperature = 0.0;  // K
    double magneticField = 0.0;       // polar field in gauss (neutron stars, white dwarfs)
    double magneticTiltDeg = 30.0;    // angle between magnetic and spin axes
    bool pulsar = false;              // beams along the magnetic axis
    bool magnetar = false;            // surface flares

    bool hasHorizon() const { return kind == Kind::BlackHole; }
    // The spin parameter a = J / M of the Kerr metric used outside: chi itself for a
    // black hole (up to 0.999); for a body capped at 0.99, since a Kerr metric with
    // a > 1 is a naked singularity, not the exterior of a planet (whose frame
    // dragging, ~2J/r^3, is negligible anyway).
    double metricSpin() const;
    // Where things stop: the event horizon r+, or the surface R.
    double boundaryRadius() const;
    // Angular velocity of the surface (d phi / dt): J / I for a body, the horizon's
    // Omega_H = a / (2 r+) for a black hole.
    double surfaceOmega() const;
    // Fastest rigid rotation before the equator sheds mass, as a spin chi.
    double maxChi() const;
    // Compactness M / R (0 for a black hole).
    double compactness() const { return hasHorizon() ? 0.0 : 1.0 / radius; }
    // Lengths the scene is laid out in, in M: 1 for compact objects, R / 6 for big
    // ones, so that every scene looks like a neutron-star-sized one on screen and the
    // GPU works with numbers of order 1-1000 (see README, "Scene scale").
    double sceneScale() const;
};

// The radius below which every object of this kind is drawn, in M: the boundary
// radius, or r+ for a black hole. (Convenience for the many call sites that only
// know a spin and an optional surface radius.)
double boundaryRadius(double a, double surface);

// --- Catalogue ------------------------------------------------------------------

struct CatalogEntry {
    const char* label;
    const char* description;
    Object object;
    double cameraDistance;  // in scene units (sceneScale() M)
    bool disk;              // shows an accretion disk by default
};
const std::vector<CatalogEntry>& catalog();

// --- Surface physics (Schwarzschild exterior unless noted) ------------------------

double surfaceRedshift(double R);         // z: 1 + z = (1 - 2/R)^(-1/2)
double surfaceGravity(double R);          // proper acceleration of a static observer, 1/M
double escapeSpeed(double R);             // locally measured, in c
double surfaceClockRate(double R);        // d tau / dt of a static clock
bool photonSphereVisible(double R);       // R < 3: the photon ring forms outside the surface
bool iscoOutside(double a, double R);     // the ISCO lies above the surface

// A rigidly rotating surface element at (R, theta) with angular velocity Omega
// (Kerr exterior): u = u^t (1, 0, 0, Omega), u^t = [-(g_tt + 2 Omega g_tphi + Omega^2 g_phiphi)]^(-1/2).
// Returns 0 if no such timelike motion exists.
double rotatingSurfaceUt(double a, double R, double theta, double omega);

// Light leaving the surface at angle alpha to the radial direction (static frame)
// reaches infinity having swept psi around the star (exact integral).
double bendingAngle(double R, double alpha);
// Beloborodov (2002): 1 - cos(alpha) = (1 - cos(psi)) (1 - 2/R), good to ~1% for R > 4.
double bendingAngleBeloborodov(double R, double alpha);
// The angle reached by light leaving tangentially: points of the surface farther
// than this from the line of sight are hidden. pi or more: the whole surface shows.
double maxBendingAngle(double R);
double visibleFraction(double R);   // fraction of the surface area seen at once
double fullVisibilityRadius();      // below this R the whole surface is visible (~3.52)

// Exact Schwarzschild periapsis advance per orbit for semi-latus rectum p (in M) and
// eccentricity e: 4 sqrt(p / (p - 6 + 2e)) K(4e / (p - 6 + 2e)) - 2 pi.
double periapsisAdvance(double p, double e);

// --- Accretion -----------------------------------------------------------------------

// Fraction of rest mass released by matter spiralling in through a thin disk to the
// ISCO and crossing the horizon: 1 - E_isco (prograde).
double blackHoleEfficiency(double a);
struct SurfaceAccretion {
    double rInner;        // where the disk ends: max(ISCO, R)
    double disk;          // released in the disk: 1 - E(rInner)
    double boundaryLayer; // released settling onto the surface: E(rInner) - E_surface
    double total;         // 1 - E_surface
};
// Matter ends at rest on a surface rotating with Omega (co-rotating at the equator).
SurfaceAccretion surfaceEfficiency(double a, double R, double omega);

// Alfven radius where a dipole field (surface field B, gauss) stops a disk feeding
// mdot (g/s): r_A = (mu^4 / (2 G M mdot^2))^(1/7) with dipole moment mu = B R^3.
// Returns r_A in km. The disk is truncated near half of it.
double alfvenRadiusKm(double fieldGauss, double radiusKm, double massSolar, double mdotGramsPerSecond);
// Eddington accretion rate, L_Edd / (0.1 c^2), in g/s.
double eddingtonRate(double massSolar);

// --- Tidal disruption -----------------------------------------------------------------

// Where a star of radius rStar (in M) and mass ratio m = m* / M is torn apart:
// r_t = rStar m^(-1/3), in M.
double tidalRadius(double rStar, double massRatio);
// The heaviest black hole that disrupts a star of radius rStarKm and mass mStarSolar
// outside its horizon (Hills mass): (R* c^2 / 2G)^(3/2) / m*^(1/2), in M_sun.
double hillsMass(double rStarKm, double mStarSolar);

// --- Stellar structure ------------------------------------------------------------------

// The Tolman-Oppenheimer-Volkoff equations for a polytrope P = K rho0^Gamma,
// epsilon = rho0 + P / (Gamma - 1), in G = c = M_sun = 1 (lengths in 1.4766 km).
struct TovStar {
    double mass;    // gravitational mass, M_sun
    double radius;  // circumferential radius, G M_sun / c^2
};
TovStar solveTov(double K, double gamma, double centralDensity);
// The most massive stable star of the polytrope, and the central density giving it.
TovStar maxTovMass(double K, double gamma, double* centralDensity = nullptr);
// The neutron-star equation of state used here: K = 100, Gamma = 2 (the standard
// Einstein Toolkit test star; maximum mass ~1.64 M_sun).
constexpr double kNeutronStarK = 100.0, kNeutronStarGamma = 2.0;
double neutronStarMaxMass();
// Radius (km) on the stable branch for a given mass, or 0 above the maximum mass.
double neutronStarRadiusKm(double massSolar);

// White dwarf mass-radius relation (Nauenberg 1972, mu_e = 2):
// R = 0.0112 R_sun sqrt((M / M_Ch)^(-2/3) - (M / M_Ch)^(2/3)), in km (0 above M_Ch).
double whiteDwarfRadiusKm(double massSolar);

// --- Collapse ------------------------------------------------------------------------------

// Checks whether the object can still hold itself up; if not, turns it into what it
// collapses to (keeping M and J) and describes what happened. A star smaller than
// the Buchdahl limit, or a neutron star above the maximum mass, becomes a black hole.
// A white dwarf pushed over the Chandrasekhar mass collapses to a neutron star
// (accretion-induced collapse), or to a black hole if it is too heavy for one.
// A black hole cannot spin faster than a = 1: the excess angular momentum is shed
// first (spin capped at 0.998, Thorne's limit). Returns false if nothing happened.
bool collapseIfUnstable(Object& object, std::string& message);

// Applies the mass-shedding limit to the spin: returns true (and caps chi) if the
// requested spin was faster than the surface can turn.
bool capSpin(Object& object);

// What a star of this compactness most resembles, for the compactness slider:
// "main-sequence star", "white dwarf", "neutron star".
Kind stageForCompactness(double R, Kind original);
// Turns the object into a star of that kind (keeping M, R and J), with that kind's
// typical surface temperature and field. Returns false if it already was one.
bool becomeStage(Object& object, Kind stage, double originalTemperature);

// Neutron stars are ~1e6 K and their polar caps ~1e7 K: they shine in X-rays. They are
// drawn with a compressed color temperature, T for T <= 30,000 K and
// 30,000 (T / 30,000)^0.3 above, so that hotter still looks bluer.
double displayTemperature(double kelvin);

}  // namespace central
