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

#include <array>
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
constexpr double kSunTemperature = 5772.0;            // K
constexpr double kEarthRadiusKm = 6371.0;
constexpr double kEarthMassSolar = 3.003e-6;
constexpr double kJupiterRadiusKm = 69911.0;
constexpr double kMoonRadiusKm = 1737.4;

double lengthKm(double lengthM, double massSolar);
// A length in km as a multiple of G M / c^2.
double kmToM(double km, double massSolar);
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
    // A young neutron star drawn as films show one: a blazing core in its glow, two
    // needle-thin beams sweeping round, and starquakes that set its crust ringing.
    bool beacon = false;

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

// What sort of place the center makes, which decides what can really be found around
// it (a Moon around the Earth, planets around the Sun, plasma around a neutron star).
enum class Host { BlackHole, NeutronStar, WhiteDwarf, SunLike, RedDwarf, BrownDwarf, Giant, Earthlike, GasGiant };
Host hostOf(const Object& object);
constexpr unsigned hostBit(Host host) { return 1u << static_cast<unsigned>(host); }

// Spin chi = c J / (G M^2) = I Omega / M^2 of a body turning once every periodSeconds,
// with moment of inertia inertia * M R^2.
double chiFromPeriod(double periodSeconds, double massSolar, double radiusKm, double inertia);

// --- Catalogue ------------------------------------------------------------------

struct CatalogEntry {
    const char* label;
    const char* description;
    Object object;
    double cameraDistance;  // in scene units (sceneScale() M)
    bool disk;              // shows an accretion disk by default
    double diskOuter = 0.0; // its outer edge in surface radii (0: the scene default)
    double cameraInclinationDeg = 0.0;  // the view's angle from the spin axis (0: the default)
    double timeScale = 0.0; // pace of the simulation, scene time per second (0: the default)
    bool jets = false;      // shows jets by default (X-ray binaries launch them)
    double cameraRollDeg = 0.0;  // the picture turned about the view (0: the spin axis upright)
};
const std::vector<CatalogEntry>& catalog();

// --- Sizes ---------------------------------------------------------------------------

// A familiar thing of about the same size, for the readouts: "about the size of a city
// (Manhattan is 22 km long)", "1.4 times the Earth", "100 times the Sun".
std::string sizeComparison(double radiusKm);

// --- Light -----------------------------------------------------------------------------

// Sunlight falling on a body at distanceKm from a star of radius radiusKm and surface
// temperature T, relative to the Sun's at the Earth: (R/R_sun)^2 (T/T_sun)^4 (1 AU / d)^2.
double fluxRelativeToEarth(double radiusKm, double temperature, double distanceKm);

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
// Hottest point of a thin disk feeding mdot (g/s) onto a mass massSolar through an inner
// edge rInnerKm (Shakura and Sunyaev): sigma T^4 = 3 G M mdot / (8 pi r^3) (1 - sqrt(r_in / r)),
// highest at r = (49/36) r_in, where T = 0.488 [3 G M mdot / (8 pi sigma r_in^3)]^(1/4). In K.
// About 1e7 K for a neutron star (X-rays), a few 1e4 K for a white dwarf (ultraviolet).
double diskPeakTemperature(double massSolar, double mdotGramsPerSecond, double rInnerKm);

// --- Pulsars -----------------------------------------------------------------------------

// The light cylinder, where co-rotation would reach the speed of light: c P / (2 pi), km.
// Field lines reaching past it are dragged open; the rest close back on the star.
double lightCylinderKm(double periodSeconds);
// Height of a pulsar's radio emission above the star, from the widths of pulse profiles
// (Kijak and Gil 2003, at 1 GHz): ~400 km P^0.3 (P in seconds), at most the light cylinder.
double radioEmissionHeightKm(double periodSeconds);
// Half-opening angle of the beam (radians): the direction of the last open field line
// where the radio is made. A dipole line at colatitude theta points theta + atan(tan(theta)/2)
// from the axis (3/2 theta near the pole), and the last open line has
// sin^2(theta) = r_em / R_LC; at most 90 deg. About 5.75 deg / sqrt(P / 1 s), as measured
// (Rankin 1993).
double beamHalfAngle(double periodSeconds);
// Half-width at half maximum of a pulsar's core beam, the bright needle at its center
// (radians): 1.22 deg / sqrt(P / 1 s) (Rankin 1990), a fifth of the cone's.
double beamCoreHalfWidth(double periodSeconds);
// Torsional oscillations of a neutron star's crust (shear waves, no radial motion), set
// ringing by a starquake: fundamental frequency of degree l, in Hz,
// f_l = f_2 sqrt((l - 1)(l + 2) / 4) with f_2 = 22 Hz for a 1.4 solar-mass star
// (Hansen and Cioffi 1980; Samuelsson and Andersson 2007), as seen in the 18-30 Hz
// oscillations after SGR 1806-20's giant flare (Israel et al. 2005).
double crustTorsionalFrequency(int l);

// --- Close binaries ------------------------------------------------------------------------

// Eggleton's (1983) Roche-lobe radius (of a sphere of the lobe's volume) over the
// separation, for the star whose mass is q times its companion's:
// 0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3))), within 1% for any q.
double rocheLobeFraction(double q);
// The inner Lagrange point L1: its distance from the star of mass M1, over the
// separation, for q = M2 / M1. The root between the stars of the Roche potential's
// gradient, -1/x^2 + q/(1 - x)^2 + (1 + q) x - q = 0 (units a = 1, G (M1 + M2) = Omega^2).
double innerLagrangeFraction(double q);
// The largest disk the companion's tide lets the star of mass M1 keep, over the
// separation: 0.6 / (1 + q) (Paczynski 1977, for 0.03 < q < 1).
double tidalTruncationFraction(double q);
// The gas stream from L1 toward the star of mass M1, ballistic in the rotating Roche
// potential (Lubow and Shu 1975), leaving L1 from rest. Points (x, y) in units of the
// separation, centered on M1, the companion on +x and the orbit turning counter-clockwise
// (Coriolis bends the stream toward +y). It stops on reaching radius `stop` from M1, or at
// its closest approach; `closest` returns that approach. At most `count` points, evenly
// spaced in time.
std::vector<std::array<double, 2>> rocheStream(double q, double stop, int count, double& closest);

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
