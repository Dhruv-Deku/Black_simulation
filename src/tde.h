#pragma once

// What follows a tidal disruption (BodySystem::disrupt tears the star apart):
//
//   fallback  The debris leaves pericenter with a spread of orbital energies. The
//             unbound part flies away; each bound parcel swings out and returns to
//             pericenter one radial period later. The rate at which mass comes back
//             rises steeply once the most bound debris returns, peaks, then fades as
//             t^(-5/3) (Rees 1988; Phinney 1989; Evans & Kochanek 1989).
//   flare     The returning stream collides with itself (relativistic apsidal
//             precession makes it cross), shocks and settles into a small disk at
//             about twice the pericenter, which accretes. For a Sun-like star and a
//             hole of a million Suns this would release ~100 times the Eddington
//             luminosity; most of that light is trapped in the inflow and what escapes
//             grows only as the logarithm (photon trapping; Begelman 1979).
//   outflow   A large share of the returning gas is blown out instead, by the shocks
//             and the trapped radiation: a wind at ~0.01-0.1 c that is wide while
//             the fallback is strong and narrows into a funnel along the axis as it
//             fades (Strubbe & Quataert 2009; Metzger & Stone 2016; Dai et al. 2018;
//             Lu & Bonnerot 2020). Seen in AT2019qiz (Nicholl et al. 2020), the event
//             behind ESO's animation of a star being spaghettified.
//
// Units: G = c = 1, lengths and times in M (the hole's mass), masses in M, unless
// stated otherwise.

#include "debris.h"
#include "kerr.h"

#include <array>
#include <vector>

namespace scene {
struct Access;  // saves and restores what is private here (src/scene_file.cpp)
}

namespace tde {

// --- Order-of-magnitude scales (Newtonian, any units with G = 1) ------------------

// A star of mass m and radius rStar torn apart by a hole of mass M. The debris keeps
// the spread of specific orbital energy across the star at the tidal radius ("frozen
// in" there; Stone, Sari & Loeb 2013), dE = M rStar / r_t^2.
struct Scales {
    double tidalRadius = 0.0;   // r_t = rStar (M / m)^(1/3)
    double energySpread = 0.0;  // dE
    double aMin = 0.0;          // semi-major axis of the most bound debris, M / (2 dE)
    double tMin = 0.0;          // its period, 2 pi sqrt(aMin^3 / M): when the first debris returns
    double peakRate = 0.0;      // m / (3 tMin): the peak of dM/dt = (m / 3 tMin) (t / tMin)^(-5/3)
};
Scales classicalScales(double M, double m, double rStar);

// The same for a star of mStar Suns and rStarKm, torn apart by a hole of holeSolar Suns,
// in physical units: tMin in seconds, peakRate in solar masses per year.
Scales physicalScales(double holeSolar, double mStar, double rStarKm);

// Period of a Kepler orbit of specific orbital energy e < 0 around mass M (G = 1):
// 2 pi M / (-2 e)^(3/2).
double keplerPeriod(double M, double e);

// --- Radial motion of a bound Kerr geodesic (M = 1, mu = 1) -------------------------

// Its turning points rp < ra (zeros of R(r)) and its radial period: the coordinate
// time between two pericenter passages, 2 integral dt/dr dr from rp to ra. With
// r = (ra + rp)/2 - (ra - rp)/2 cos(chi) the square-root singularities of dt/dr at
// the turning points cancel against dr, so a plain quadrature in chi is accurate.
// The polar motion enters through the average of sin^2(theta) (exact for Q = 0).
struct RadialOrbit {
    enum Kind { Bound, Plunges, Unbound } kind = Bound;
    double rp = 0.0, ra = 0.0;
    double period = 0.0;        // coordinate time from pericenter to pericenter
    double toPericenter = 0.0;  // Kerr-Schild time until the next pericenter passage
};
RadialOrbit radialOrbit(double a, const kerr::Constants& c, const kerr::State& s);

// --- Fallback --------------------------------------------------------------------------

// The rate at which a disruption's debris returns to pericenter, built from the
// parcels themselves: each returns after its own radial period. Rates are fractions
// of the star per unit time.
class Fallback {
public:
    // returnTimes (after the disruption) and weights of the parcels that come back;
    // totalWeight counts every parcel, unbound ones included.
    void build(const std::vector<double>& returnTimes, const std::vector<double>& weights, double totalWeight);
    double rate(double t) const;      // at t after the disruption
    double returned(double t) const;  // fraction of the star returned by then
    bool empty() const { return m_rate.empty(); }

    double boundFraction = 0.0;
    double firstReturn = 0.0;  // the most bound debris
    double peakTime = 0.0;
    double peakRate = 0.0;
    double lateSlope = 0.0;    // d ln(rate) / d ln(t) well after the peak (-5/3 for a flat dM/dE)

private:
    friend struct scene::Access;
    // Tabulated on a uniform grid in ln t, smoothed over a few percent in time.
    double m_logT0 = 0.0, m_dLogT = 1.0;
    std::vector<double> m_rate, m_cumulative;
};

// --- Flare and outflow -------------------------------------------------------------------

// Eddington luminosity of a hole of massSolar Suns (electron scattering, hydrogen), erg/s.
double eddingtonLuminosity(double massSolar);
// Luminosity that escapes, in Eddington units, when accretion releases lacc Eddington
// luminosities: all of it below Eddington, 1 + ln(lacc) above (photon trapping).
double escapingLuminosity(double lacc);

// A radially coasting wind of terminal speed vInf around mass M: speed at R from
// energy conservation, v^2 = vInf^2 + 2M/R (Newtonian; the wind starts well outside
// the hole), and the time it takes from R0 to R, the integral of dR / v(R), in closed
// form: with k = vInf^2, T(R) = sqrt(R (kR + 2M)) / k - (2M / k^1.5) ln(sqrt(kR) + sqrt(kR + 2M)).
double windSpeed(double M, double vInf, double R);
double windTravelTime(double M, double vInf, double R0, double R);
// How far the wind's leading edge has got t after leaving R0 (inverts windTravelTime).
double windFront(double M, double vInf, double R0, double t);

// Opening half-angle (radians) of the outflow launched once a fraction `formed` of the
// bound debris has come back. At first no disk has formed and the stream's self-crossing
// shock throws gas out over a wide range of angles (`first`); as a disk builds up, it
// and the radiation it traps collimate the wind into a funnel (`late`). The change
// e-folds once a fifth of the bound debris is back.
double openingAngle(double late, double first, double formed);

// One disruption, as analysed when it happens.
struct Event {
    double time = 0.0;                     // coordinate time of the disruption
    std::array<double, 3> axis{0, 0, 1};   // the debris's orbital angular momentum (unit, Cartesian)
    double pericenter = 0.0;               // of the star's orbit
    double circularization = 0.0;          // where the returning gas settles: about 2 r_p
    double starMass = 0.0;                 // in M
    double starRadius = 0.0;               // in M
    double tidalRadius = 0.0;              // in M
    double unboundFraction = 0.0;
    double plungingFraction = 0.0;         // falls straight into the hole
    Fallback fallback;
};

// Analyses the debris of a star (its parcels' constants and positions) torn apart at
// `time` around a hole of spin a. At most maxSamples parcels are used (evenly spread).
Event analyze(double a, double time, const std::vector<ParticleInit>& debris, int maxSamples = 30000);

// How the outflow is launched.
struct OutflowParams {
    double fraction = 0.5;      // share of the returning debris blown out; the rest is accreted
    double speed = 0.1;         // terminal speed on the axis (c)
    double openingLate = 0.38;  // half-opening angle once the fallback has faded (radians)
    double widening = 2.4;      // ... times this at first, before a disk has formed
};

// Everything the disruptions so far add up to, at any coordinate time t. Rates add up
// over the events; the geometry (axis, launch radius) is the strongest event's.
class Aftermath {
public:
    // holeSolar sets the physical scale of the flare: what a Sun-like star torn apart by
    // this hole would release (the simulated star is scaled up to be seen).
    Aftermath(const std::vector<Event>& events, const OutflowParams& params, double holeSolar, double spin);

    bool active() const { return m_reference >= 0; }
    const Event& reference() const { return m_events[static_cast<size_t>(m_reference)]; }
    const OutflowParams& params() const { return m_params; }

    double level(double t) const;     // fallback rate relative to its peak
    double returned(double t) const;  // star masses (of the strongest event's star) returned by t
    double outflowRate(double t) const { return m_params.fraction * level(t); }
    double opening(double t) const;   // half-opening angle of the outflow launched at t
    double luminosity(double t) const;       // escaping luminosity relative to its peak
    double eddingtonRatio(double t) const;   // ... in Eddington units (Sun-like star)
    double accretionRatio(double t) const;   // accretion power in Eddington units, before trapping
    double start() const { return m_start; }  // when the returning debris starts to arrive (1% of the peak)
    double launchRadius() const;             // the circularization radius, at least 10 M
    double front(double t) const;            // radius of the outflow's leading edge (launch radius before)
    double peakTime() const { return m_peakTime; }

private:
    std::vector<Event> m_events;
    int m_reference = -1;
    OutflowParams m_params;
    double m_peakMass = 0.0;   // peak of the summed rate, star masses per M
    double m_peakAccretion = 0.0;  // accretion power at the peak, Eddington units
    double m_start = 0.0, m_peakTime = 0.0;
};

}  // namespace tde
