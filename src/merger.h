#pragma once

// Two neutron stars spiral together, merge and collapse into a black hole (src/merger_sim.*
// runs the stars as smoothed-particle hydrodynamics on the GPU; this file holds the physics
// that does not need a GPU, and is tested in tests/test_kerr.cpp).
//
//   inspiral   The binary loses energy to gravitational waves and its orbit shrinks at the
//              rate Peters (1964) derived: a(t) = a0 (1 - t / T)^(1/4), T = 5 a0^4 / (256 m1 m2 M).
//              The waves come at twice the orbital frequency and sweep up in a "chirp" whose
//              shape is set by the chirp mass alone (GW170817: 1.188 Suns, 30 Hz to ~1.5 kHz in
//              the last 100 s).
//   merger     Tides stretch the stars just before they touch; their cores plunge together,
//              shock-heat to tens of MeV and spin into a hypermassive neutron star, held up
//              for a while by its rotation and heat, flinging out spiral arms of neutron-rich
//              matter (Shibata & Uryu 2000; Rosswog et al. 1999; Bauswein & Janka 2012).
//   collapse   The remnant collapses into a black hole of ~2.6 Suns spinning at ~0.7, which
//              rings down at its quasi-normal frequency (Berti, Cardoso & Will 2006); the
//              gas with too much angular momentum to fall in stays behind as a torus.
//   GRB        A jet launched along the spin axis drills through the ejecta and breaks out:
//              a short gamma-ray burst (GRB 170817A came 1.7 s after the merger).
//   kilonova   The ejecta (a few hundredths of a Sun at 0.1-0.3 c) are heated by the decay of
//              freshly made r-process nuclei (Li & Paczynski 1998; Metzger et al. 2010) and
//              shine for days: a blue, lanthanide-poor component first, then a red one
//              (Villar et al. 2017's three-component fit of AT2017gfo).
//
// Units: G = c = 1 and the Sun's mass is 1, so lengths are in G M_sun / c^2 = 1.476625 km
// and times in G M_sun / c^3 = 4.925491 microseconds, unless stated otherwise.

#include <glm/glm.hpp>

#include <array>
#include <vector>

namespace merger {

constexpr double kPi = 3.14159265358979323846;
constexpr double kKm = 1.476625;            // km per length unit
constexpr double kSecond = 4.925491e-6;     // seconds per time unit
constexpr double kGramsPerCc = 6.17714e17;  // g/cm^3 per density unit, M_sun / (G M_sun / c^2)^3
constexpr double kSunGrams = 1.98847e33;
constexpr double kLightCm = 2.99792458e10;  // cm/s
constexpr double kDay = 86400.0;            // s

// --- Inspiral -----------------------------------------------------------------------------

double chirpMass(double m1, double m2);  // (m1 m2)^(3/5) / (m1 + m2)^(1/5)
// Time for a circular binary of separation a to coalesce (Peters 1964), 5 a^4 / (256 m1 m2 M).
double coalescenceTime(double a, double m1, double m2);
// Its separation a time t later, a (1 - t / T)^(1/4) (0 once it has coalesced).
double separationAfter(double a, double t, double m1, double m2);
// Gravitational-wave frequency (twice the orbital one, cycles per unit time) a time tau
// before coalescence: (1/pi) (5 / (256 tau))^(3/8) Mc^(-5/8).
double chirpFrequency(double tau, double mc);
// Power radiated by a circular binary, (32/5) eta^2 (M / a)^5 (units of c^5/G).
double gwLuminosity(double a, double m1, double m2);

// The 2.5PN radiation-reaction acceleration of the relative orbit (x = x1 - x2, v = v1 - v2):
// (8/5) eta M^2 / r^3 [ (3 v^2 + 17/3 M/r) rdot n - (v^2 + 3 M/r) v ] (Iyer & Will 1993, in
// harmonic coordinates). On a circular orbit it drains energy at exactly Peters' rate.
glm::dvec3 radiationReaction(const glm::dvec3& x, const glm::dvec3& v, double m1, double m2);

// The wave seen from direction n (unit, from the source toward the observer) at distance D:
// the transverse-traceless part of h_ij = (2 / D) d^2 I_ij / dt^2 (quadrupole formula),
// on the polarization axes p = e_theta, q = e_phi of n. iddot holds d^2 I / dt^2 as
// xx, yy, zz, xy, xz, yz.
struct Strain {
    double plus = 0.0, cross = 0.0;
};
Strain strain(const std::array<double, 6>& iddot, const glm::dvec3& n, double distance);
// d^2 I / dt^2 of point masses: sum of m (2 v_i v_j + x_i a_j + x_j a_i).
std::array<double, 6> quadrupoleSecondDerivative(const std::vector<glm::dvec3>& x, const std::vector<glm::dvec3>& v,
                                                 const std::vector<glm::dvec3>& a, const std::vector<double>& m);

// --- The stars: n = 1 polytropes ------------------------------------------------------------

// P = K rho^2. The radius of an n = 1 polytrope does not depend on its mass,
// R = sqrt(pi K / 2) (G = 1): K = 2 R^2 / pi.
double polytropeK(double radius);
double polytropeCentralDensity(double mass, double radius);  // pi M / (4 R^3)
double polytropeDensity(double mass, double radius, double r);  // rho_c sin(xi) / xi, xi = pi r / R
double polytropeEnclosed(double mass, double radius, double r);  // M (sin xi - xi cos xi) / pi
// The radius inside which a fraction f of the mass lies.
double polytropeRadiusOfFraction(double radius, double f);
// n particle positions filling a star of this radius with equal masses: a close-packed
// lattice inside the unit sphere, each point moved outward so that the share of points
// inside it equals the share of the polytrope's mass (a "stretched lattice": far quieter
// than random positions, so the star starts close to equilibrium).
std::vector<glm::dvec3> starParticles(int n, double radius);

// Thermal part of the hybrid equation of state the merger runs with:
// P = K rho^2 + (gammaThermal - 1) rho u_th, gammaThermal = 1.75 (Bauswein et al. 2010).
constexpr double kGammaThermal = 1.75;
// Temperature of nuclear matter with thermal energy u_th per unit rest mass (an ideal
// gas of nucleons: u_th = (3/2) kT / m_n c^2), in kelvin. A rough guide: degenerate
// matter is hotter for the same u_th, and pairs and photons share the heat.
double temperatureFromThermalEnergy(double uThermal);
double thermalEnergyFromTemperature(double kelvin);

// --- Collapse ---------------------------------------------------------------------------------

// The black hole the remnant collapses into. Walking outward from the densest point, a
// bound parcel falls in unless its specific angular momentum exceeds that of the innermost
// stable circular orbit of the mass inside it, 2 sqrt(3) M (Schwarzschild). The rest stays
// behind: a torus, and what is already unbound.
struct BlackHole {
    double mass = 0.0;
    double spin = 0.0;           // chi = J / M^2 (capped at 0.998, Thorne's limit)
    glm::dvec3 position{0.0};
    glm::dvec3 velocity{0.0};
    glm::dvec3 axis{0, 0, 1};    // its spin axis
    double diskMass = 0.0;       // bound mass left outside
    double unboundMass = 0.0;
};
// x, v, m of the particles, their specific potential energy phi (to tell bound from
// unbound) and which may be captured (alive); swallowed is set for those that fall in.
BlackHole collapse(const std::vector<glm::dvec3>& x, const std::vector<glm::dvec3>& v, const std::vector<double>& m,
                   const std::vector<double>& phi, const std::vector<bool>& alive, const glm::dvec3& center,
                   std::vector<bool>& swallowed);

// The ringing of the new hole: the l = m = 2 fundamental quasi-normal mode (Berti, Cardoso
// & Will 2006 fits): M omega = 1.5251 - 1.1568 (1 - chi)^0.1292, quality factor
// Q = 0.7 + 1.4187 (1 - chi)^-0.4990, damping time 2 Q / omega.
struct Ringdown {
    double omega = 0.0;   // angular frequency (rad per unit time)
    double tau = 0.0;     // e-folding time of the amplitude
};
Ringdown ringdown(double mass, double spin);

// Kerr horizon radius and Schwarzschild-like shadow (critical impact parameter, sqrt(27) M).
double horizonRadius(double mass, double spin);

// --- The debris cloud -------------------------------------------------------------------------

// What the collision itself throws out, the dynamical ejecta (Hotokezaka et al. 2013;
// Sekiguchi et al. 2016; Radice et al. 2018): tidal tails flung from the stars' far sides
// into the orbital plane, and gas squeezed out of the shocked interface between the cores,
// faster and more toward the poles. A few thousandths of a Sun at 0.1-0.3 c, launched within
// about a millisecond of contact: a cloud bursts out of the remnant as the stars touch,
// coasts outward (r = v t) and cools as it expands. Later the remnant's and torus's winds
// add most of the kilonova's mass (the components below).
struct DebrisCloud {
    double mass = 0.005;               // solar masses (GW170817's fits: ~0.001-0.01 dynamical)
    double edgeSpeed = 0.3;            // c: its fastest gas, the shock-driven part
    double launchKm = 15.0;            // where it starts from: the merging stars' surface
    double launchTemperature = 5e10;   // K: shocked to a few MeV where the cores collide
    double opacity = 0.08;             // cm^2/g: electron scattering, ~0.4 Ye with Ye ~ 0.2
};
// Its edge t seconds after contact, km: launchKm + edgeSpeed c t.
double debrisRadiusKm(const DebrisCloud& c, double t);
// Its temperature then, K. The hot ejecta are dominated by radiation and pairs and expand
// adiabatically, T^3 V constant, so T falls as 1 / R.
double debrisTemperature(const DebrisCloud& c, double t);
// Its optical depth from the centre to the edge, 3 kappa M / (4 pi R^2) for a uniform
// sphere: ~10^15 at first, and still hundreds after an hour to electron scattering alone
// (once it cools, the r-process elements' lines keep it opaque for days: the kilonova).
double debrisOpticalDepth(const DebrisCloud& c, double t);

// --- Kilonova ---------------------------------------------------------------------------------

// One component of the ejecta, coasting homologously (r = v t).
struct EjectaComponent {
    const char* name = "";
    double mass = 0.01;         // solar masses
    double velocity = 0.2;      // c (mean)
    double opacity = 1.0;       // cm^2/g: low without lanthanides, ~10 with
    double temperatureFloor = 1000.0;  // K: the photosphere recedes once it cools to this
};
// Villar et al. (2017), the three-component fit of GW170817's kilonova: blue (lanthanide-
// poor, fast, polar), purple, red (lanthanide-rich).
std::array<EjectaComponent, 3> gw170817Ejecta();

// Specific heating rate of r-process matter (Korobkin et al. 2012), erg/g/s, t in s:
// 4e18 (0.5 - arctan((t - 1.3) / 0.11) / pi)^1.3.
double rProcessHeating(double t);
// The share of that heat the ejecta keep (Barnes et al. 2016's fit): the rest escapes as
// neutrinos and as gamma-rays the thinning ejecta no longer stop. t in s.
double thermalization(double t);

// Light curve of one component (Arnett's model as Villar et al. 2017 use it): luminosity
// L(t) = exp(-(t/td)^2) integral_0^t 2 M eps(t') (t'/td) exp((t'/td)^2) dt' / td, with
// diffusion time td = sqrt(2 kappa M / (beta v c)), beta = 13.7; photosphere temperature
// (L / (4 pi sigma (v t)^2))^(1/4) until it cools to the floor, then held there while the
// photosphere recedes. Tabulated on a log grid in t.
class Kilonova {
public:
    explicit Kilonova(const EjectaComponent& c);
    double luminosity(double t) const;   // erg/s, t in s since the merger
    double temperature(double t) const;  // K
    double photosphere(double t) const;  // cm
    double peakTime() const { return m_peakTime; }
    double peakLuminosity() const { return m_peakLuminosity; }
    double diffusionTime() const { return m_td; }
    const EjectaComponent& component() const { return m_c; }

private:
    double sample(const std::vector<double>& table, double t) const;
    EjectaComponent m_c;
    double m_td = 1.0;
    double m_logT0 = -3.0, m_dLogT = 0.01;
    std::vector<double> m_L, m_T, m_R;
    double m_peakTime = 0.0, m_peakLuminosity = 0.0;
};

// --- The jet (a short gamma-ray burst) ----------------------------------------------------------

// A jet launched along the spin axis tlaunch after the merger drills through the ejecta
// with its head at headSpeed (c). The ejecta's edge coasts at vmax (the blue component's
// along the axis, kEjectaEdge x kBluePolarStretch x 0.256 c = 0.37 c), so the head catches it at
// tbreakout = headSpeed tlaunch / (headSpeed - vmax); the gamma-rays escape then.
// The ejecta's edge: each component reaches out to this times its mean speed, stretched
// into its shape (shaders/merger_view.frag, ejectaStretch): the blue toward the poles (x 1.15
// along the axis, 0.8 in the plane), the red into a flattened disk (x 1.6 in the plane,
// 0.55 along the axis), the purple slightly flattened (1.15, 0.85).
constexpr double kEjectaEdge = 1.25;
constexpr double kBluePolarStretch = 1.15;

struct Jet {
    double launch = 0.8;       // s after the merger
    double headSpeed = 0.7;    // c, while it drills through the ejecta
    double coreAngle = 0.09;   // radians: the narrow, ultra-relativistic core (~5 degrees)
    double lorentz = 100.0;    // core Lorentz factor
    double cocoonSpeed = 0.3;  // c: the hot cocoon spills sideways after breakout
};
double breakoutTime(const Jet& jet, double vmax);
// Doppler factor 1 / (Gamma (1 - beta cos angle)) of the jet's core seen at angle from
// its axis: why a jet seen off-axis looks faint (GRB 170817A was seen ~20-30 degrees off).
double dopplerFactor(double lorentz, double angle);

// --- How time is shown ------------------------------------------------------------------------

// The event spans milliseconds (the merger) to weeks (the kilonova). Time is shown at a
// rate (physical seconds per second of viewing) that grows with the time to or from the
// merger, so each decade takes about as long to watch, and slows down for the merger, the
// collapse and the burst.
struct Tempo {
    double inspiralRate = 0.32;  // e-folds of the time to merger per second of viewing
    double mergerRate = 0.0004;  // s per s at the merger itself (2,500 times slow motion)
    double afterRate = 0.18;     // e-folds per second of viewing after the merger
    double lateRate = 0.55;      // ... after the black hole forms
    double burstHold = 1.5;      // seconds of viewing around the gamma-ray burst
};
// Physical seconds per second of viewing at time t (s, from the merger; negative before),
// given when the hole forms and when the burst breaks out.
double tempo(const Tempo& p, double t, double collapse, double burst);

}  // namespace merger
