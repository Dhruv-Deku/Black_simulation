#pragma once

// Kerr black hole physics.
//
// Units are geometrized with G = c = M = 1, so distances are in units of the
// gravitational radius GM/c^2 and times in units of GM/c^3. The spin parameter
// a = J/M lies in [0, 1); orbits are "prograde" when they co-rotate with the hole.
//
// Coordinates are ingoing Kerr-Schild (t, r, theta, phi). r and theta are the
// Boyer-Lindquist ones; t and phi are shifted along each ingoing light ray,
//     t   = t_BL   + integral (r^2 + a^2) / Delta dr - r,
//     phi = phi_BL + integral a / Delta dr,
// which removes the Boyer-Lindquist coordinate singularity at the horizon: bodies,
// gas, light and the camera all cross it in finite time and carry on inside. Far
// from the hole the two systems agree; everything observable is the same.
//
// Geodesics are integrated with Carter's separated equations in Mino time tau,
// defined by d(lambda) = Sigma d(tau). In Mino time the radial and polar motions
// decouple, and using their second-order form
//     d^2 r / d tau^2     = R'(r) / 2
//     d^2 theta / d tau^2 = Theta'(theta) / 2
// avoids the sign flips at turning points that the first-order form needs.
// The same equations drive the GPU ray tracer and the debris particles.

#include <array>

namespace kerr {

constexpr double kPi = 3.14159265358979323846;

// --- Characteristic radii --------------------------------------------------

double horizonRadius(double a);                      // outer event horizon r+
double innerHorizonRadius(double a);                 // inner (Cauchy) horizon r-
double ergosphereRadius(double a, double theta);     // static limit
double iscoRadius(double a, bool prograde);          // innermost stable circular orbit
double photonOrbitRadius(double a, bool prograde);   // circular equatorial photon orbit
double marginallyBoundRadius(double a, bool prograde);

// How deep inside the horizon things are followed: just above the inner horizon,
// where the Kerr interior stops being trustworthy (the mass-inflation instability),
// or just above the singularity of a non-spinning hole.
double interiorLimit(double a);

// --- Metric quantities -------------------------------------------------------

struct Metric {
    double sigma;  // r^2 + a^2 cos^2(theta)
    double delta;  // r^2 - 2r + a^2
    double A;      // (r^2 + a^2)^2 - a^2 delta sin^2(theta)
    double omega;  // frame-dragging angular velocity 2ar / A
};
Metric metricAt(double a, double r, double theta);

using Vec4 = std::array<double, 4>;  // components along (t, r, theta, phi)

// Covariant Kerr-Schild metric g_mu_nu, regular everywhere except at r = 0.
struct MetricTensor {
    double g[4][4];
};
MetricTensor metricTensor(double a, double r, double theta);
double dot(const MetricTensor& g, const Vec4& x, const Vec4& y);
Vec4 lower(const MetricTensor& g, const Vec4& x);

// t - t_BL and phi - phi_BL at radius r (away from the horizons). The azimuth shift
// vanishes far from the hole; the time shift grows there like 2 ln r.
double timeShift(double a, double r);
double azimuthShift(double a, double r);

// --- Circular equatorial orbits (Bardeen, Press & Teukolsky 1972) ------------

double circularOmega(double a, double r, bool prograde);      // d(phi)/dt
double circularZamoSpeed(double a, double r, bool prograde);  // speed seen by a ZAMO

struct Constants {
    double E;   // energy per unit mass (photons: arbitrary scale)
    double L;   // axial angular momentum
    double Q;   // Carter constant
    double mu;  // rest mass: 1 for particles, 0 for photons
};
Constants circularOrbit(double a, double r, bool prograde);

// --- Geodesics in Mino time --------------------------------------------------

struct State {
    double t, r, theta, phi;
    double vr;      // dr/dtau
    double vtheta;  // dtheta/dtau
};

struct Derivatives {
    double dt, dr, dtheta, dphi, dvr, dvtheta;  // d/dtau of each State field
};

Derivatives derivatives(double a, const Constants& c, const State& s);

double radialPotential(double a, const Constants& c, double r);       // R(r)  = (dr/dtau)^2
double polarPotential(double a, const Constants& c, double theta);    // Theta = (dtheta/dtau)^2

// (P + vr) / Delta, with P = E (r^2 + a^2) - a L: the part of dt/dtau and dphi/dtau
// that is singular at the horizon in Boyer-Lindquist coordinates. For inward motion
// numerator and denominator vanish together there; on the constraint surface
// vr^2 = R = P^2 - Delta K the ratio equals K / (P - vr), which is smooth through
// the horizon. K = mu^2 r^2 + (L - a E)^2 + Q.
double inflowFactor(double P, double vr, double delta, double K);

// Photons whose closest approach to the spin axis, about |L| / sqrt(Q + a^2 E^2)
// radians, would be smaller than this are sent exactly through the axis instead (L = 0).
constexpr double kAxisSnap = 3e-4;

// Launches a geodesic from (r, theta, phi) with a velocity measured by the local
// zero-angular-momentum observer (ZAMO). v holds the (r, theta, phi) components of
// the 3-velocity in the ZAMO's orthonormal frame:
//   massive (mu = 1): |v| < 1 is the physical speed;
//   photon  (mu = 0): v is the unit propagation direction, and the result is
//                     normalized to E = 1 (and snapped to the axis, see kAxisSnap).
// ZAMOs exist only outside the horizon.
void launchFromZamo(double a, double r, double theta, double phi, const double v[3],
                    double mu, Constants& c, State& s);

// Inverse of launchFromZamo for a massive particle: the 3-velocity it has in the
// local ZAMO frame (outside the horizon).
void zamoVelocity(double a, const Constants& c, const State& s, double v[3]);

// Contravariant 4-velocity (u^t, u^r, u^theta, u^phi) of a massive particle.
Vec4 fourVelocity(double a, const Constants& c, const State& s);

// 4-velocity of the ZAMO at (r, theta), outside the horizon.
Vec4 zamoFourVelocity(double a, double r, double theta);

// --- Observers -------------------------------------------------------------------

// An observer's orthonormal frame: e[0] is its 4-velocity, e[1], e[2], e[3] its unit
// spatial axes, built to follow the r, theta and phi directions. For a ZAMO this is
// the ZAMO's own frame; for a moving observer it is the frame it carries along, in
// which the sky appears aberrated. It exists inside the horizon too.
struct Frame {
    Vec4 e[4];
};
Frame observerFrame(double a, double r, double theta, const Vec4& u);

// The photon an observer sees arriving from direction `look` (a unit vector on the
// frame's spatial axes). Fills its constants, scaled so that the observer measures
// unit energy (E can then be zero or negative inside the ergosphere), and its state
// at the observer, moving forward in time toward it. shaders/trace.frag performs the
// same computation for every pixel.
void photonFromFrame(double a, double r, double theta, double phi, const Frame& f, const double look[3],
                     Constants& c, State& s);

// The same for an observer moving with 3-velocity vObs relative to the ZAMO at
// (r, theta, phi): the photon normalized to E = 1, and the frequency the observer
// measures for it.
double observedPhoton(double a, double r, double theta, double phi, const double vObs[3], const double look[3],
                      Constants& c, State& s);

// --- Integration -------------------------------------------------------------------

// One RK4 step of size dtau in Mino time.
void stepMino(double a, const Constants& c, State& s, double dtau);

// The second-order equations conserve vr^2 - R(r) exactly, including any error in
// it. R grows like r^4, so a relative rounding error far from the hole becomes an
// O(1) absolute error near it, enough to create fake turning points. Away from
// turning points this re-projects the velocities onto vr^2 = R and
// vtheta^2 = Theta; near them it leaves the second-order dynamics alone to carry
// the motion through the bounce. Photon tracing calls it after every step.
void enforceConstraints(double a, const Constants& c, State& s);

// Mino-time step size used to trace photons; shaders/trace.frag uses the same rule.
double photonStepSize(double a, const Constants& c, const State& s, double tolerance);

// One RK4 step that advances coordinate time by exactly dt.
void stepCoordinateTime(double a, const Constants& c, State& s, double dt);

// A coordinate-time step size that keeps the per-step motion small.
double suggestedTimeStep(double a, const Constants& c, const State& s, double tolerance);

// Keeps theta in [0, pi] by reflecting through the pole (same point in space).
void normalizeAngles(State& s);

// --- Equatorial bound orbits ---------------------------------------------------

// Solves for (E, L) of the equatorial orbit with periapsis rp and apoapsis ra.
// Returns false if no such bound orbit exists.
bool equatorialOrbit(double a, double rp, double ra, bool prograde, double& E, double& L);

// Smallest periapsis for which a bound equatorial orbit with apoapsis ra exists.
// Orbits just above it "zoom" out to ra and "whirl" around the hole near rp.
bool separatrixPeriapsis(double a, double ra, bool prograde, double& rp);

// --- Cartesian coordinates -------------------------------------------------------

// Kerr-Schild Cartesian coordinates, z = spin axis:
//     x + i y = (r + i a) e^{i phi} sin(theta),  z = r cos(theta).
// In them the metric is flat space plus a single null term, g = eta + f k k.
std::array<double, 3> toCartesian(double a, double r, double theta, double phi);

// Inverse of toCartesian: (r, theta, phi).
std::array<double, 3> fromCartesian(double a, double x, double y, double z);

// Cartesian velocity of a point moving with coordinate rates (dr, dtheta, dphi).
std::array<double, 3> cartesianVelocity(double a, double r, double theta, double phi, double dr, double dtheta,
                                        double dphi);

// The local orthonormal r, theta, phi directions at a point, as Cartesian vectors.
// They match the ZAMO frame's axes up to O(a / r^2) radians.
struct LocalAxes {
    std::array<double, 3> r, theta, phi;
};
LocalAxes localAxes(double a, double r, double theta, double phi);

}  // namespace kerr
