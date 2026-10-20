#pragma once

// The Penrose process: drawing energy out of a spinning black hole.
//
// Inside the ergosphere, where even light must turn with the hole, a particle can have
// negative energy as seen from far away: moving against the rotation fast enough,
// E = gamma (alpha + omega varpi v_phi) < 0 (per unit mass, from its velocity v measured by
// the local zero-angular-momentum observer, ZAMO; alpha = sqrt(Sigma Delta / A) is the
// lapse, omega = 2 a r / A the frame dragging and varpi = sqrt(A / Sigma) sin(theta) the
// circumferential radius). That needs v_phi < -alpha / (omega varpi), possible below the
// speed of light only where alpha < omega varpi: inside the ergosphere.
//
// A body falling in splits there. Its 4-momentum is shared (in the ZAMO's frame,
// m0 gamma0 (1, v0) = m1 gamma1 (1, v1) + m2 gamma2 (1, v2)), so energy and angular
// momentum are too: m0 E0 = m1 E1 + m2 E2. With E1 < 0, fragment 2 leaves with more
// than the whole had, m2 E2 - m0 E0 = -m1 E1 > 0. Fragment 1 falls in, and the hole's
// mass and angular momentum change by dM = m1 E1 < 0, dJ = m1 L1: it spins down. Its
// area never shrinks: dM - Omega_H dJ = m1 (E1 - Omega_H L1) >= 0, Omega_H = a / (2 r+).
// What can be extracted at all is M - M_irr, with M_irr = sqrt(A / 16 pi) = sqrt(M r+ / 2)
// (Christodoulou: M^2 = M_irr^2 + J^2 / (4 M_irr^2)), 29% of M for a = 1; one split at the
// horizon of a maximally spinning hole gains at most (sqrt 2 - 1) / 2 = 20.7%.
//
// Units: G = c = M = 1; masses in M, energies and angular momenta per unit rest mass.

#include <string>

namespace penrose {

// The local ZAMO's lapse alpha, frame-dragging rate omega and circumferential radius varpi.
struct Zamo {
    double alpha, omega, varpi;
};
Zamo zamoAt(double a, double r, double theta);

// A 4-momentum on the ZAMO's orthonormal axes (t, r, theta, phi).
struct Momentum {
    double p[4] = {0.0, 0.0, 0.0, 0.0};
    double mass() const;  // sqrt(p0^2 - |p|^2), or -1 if spacelike
};
Momentum momentum(double mass, const double v[3]);  // m gamma (1, v)
// Energy -p_t = alpha p^0 + omega varpi p^phi, and angular momentum p_phi = varpi p^phi (totals).
double energy(const Zamo& z, const Momentum& p);
double angularMomentum(const Zamo& z, const Momentum& p);

struct Split {
    bool valid = false;
    std::string why;               // when not
    double r = 0.0, theta = 0.0;   // where
    double m0 = 0.0, m1 = 0.0, m2 = 0.0;  // rest masses
    double E0 = 0.0, E1 = 0.0, E2 = 0.0;  // per unit rest mass
    double L0 = 0.0, L1 = 0.0, L2 = 0.0;
    double v1[3] = {0.0, 0.0, 0.0}, v2[3] = {0.0, 0.0, 0.0};  // the fragments' ZAMO-frame velocities
    double efficiency = 0.0;       // (m2 E2 - m0 E0) / (m0 E0) = -m1 E1 / (m0 E0)
    double spinBefore = 0.0;       // the hole's a = J / M then
    bool absorbed = false;         // fragment 1 has crossed the horizon
    bool applied = false;          // ... and the hole's spin has been changed by it
};

// A body of mass m0 moving at v0 (ZAMO frame) at (r, theta) throws off fragment 1 of mass
// m1 at v1; fragment 2 takes the rest of the 4-momentum. Invalid if what is left is not
// a particle (spacelike or past-directed).
Split split(double a, double r, double theta, double m0, const double v0[3], double m1, const double v1[3]);

// The split that gains the most, throwing fragment 1 off at `speed`: its direction (in
// the r-phi plane) and mass are searched, keeping only splits where fragment 1 falls into
// the hole and fragment 2 escapes to infinity. Invalid (with the reason) if none gains.
Split bestSplit(double a, double r, double theta, double m0, const double v0[3], double speed);

// The angular momentum of the prograde parabolic orbit (E = 1, in the equator) whose
// closest approach is r: the smaller root of (2 - r) L^2 - 4 a L + 2 (r^2 + a^2) = 0, from
// R(r) = 0. False if no such orbit comes from infinity that close (r below the marginally
// bound orbit's radius).
bool parabolicPeriapsis(double a, double r, double& L);
// Where a body falling in from far away on such an orbit should turn round to split: well
// inside the ergosphere (70% of the way from its edge to the horizon), but no closer than
// a parabolic orbit can come.
double splitRadius(double a);

// The hole of mass M and spin a = J / M.
double horizonArea(double M, double a);      // 8 pi M r+ = 4 pi (r+^2 + a^2)
double irreducibleMass(double M, double a);  // sqrt(A / 16 pi)
double horizonOmega(double M, double a);     // a / (2 M r+)
// Its spin after swallowing fragment 1 of a split (a' = J' / M'^2, in units of the old M).
double spinAfter(double a, const Split& s);
double massAfter(const Split& s);            // 1 + m1 E1

// One split at the horizon of a hole with a = 1, from rest at infinity: (sqrt 2 - 1) / 2.
constexpr double kMaxEfficiency = 0.20710678118654752;

}  // namespace penrose
