#pragma once

// Tidal deformation of a star: the affine model (Carter & Luminet 1983).
//
// The star is a homogeneous, adiabatic ellipsoid whose every element moves
// together: a point X of the undisturbed star (radius R) sits at x = F X, for one
// 3x3 matrix F that evolves in the star's proper time as
//     d^2F/dtau^2 = C F - (3/2) w^2 A(F) F / det F + w^2 (det F)^(1 - gamma) F^-T
// The three terms are the hole's tide (tidal tensor C), the star's own gravity
// (A holds the index symbols of the ellipsoid; a sphere has A = 2/3) and its
// pressure, balancing exactly at F = 1. w^2 = m / R^3 is the star's dynamical
// frequency squared; for a star whose tidal radius is r_t it is 1 / r_t^3 (but see
// kCondensation).
//
// Far from the hole the star stays round. Near its tidal radius it is pulled into
// a prolate egg; on a deep pass it is squeezed into a pancake at periapsis and then
// stretched without limit: it has been disrupted, and the debris keeps the stretched
// shape and the internal velocity field dF/dtau X.

#include "kerr.h"

#include <array>

namespace tidal {

using Mat3 = std::array<std::array<double, 3>, 3>;
using Vec3 = std::array<double, 3>;

Mat3 identity();
Mat3 multiply(const Mat3& a, const Mat3& b);
Mat3 transpose(const Mat3& m);
double determinant(const Mat3& m);
Mat3 inverse(const Mat3& m);

// Carlson's symmetric elliptic integral R_D(x, y, z).
double carlsonRD(double x, double y, double z);

// Eigenvalues (ascending) and unit eigenvectors (columns) of a symmetric matrix.
void symmetricEigen(const Mat3& m, Vec3& values, Mat3& vectors);

// Index symbols A_i of an ellipsoid with semi-axes (in any units) a[3]; they sum to 2.
Vec3 indexSymbols(const Vec3& a);

// Tidal tensor of the hole felt by a body on the geodesic (c, s): a point displaced
// by x from its center accelerates by C x (Kerr-Schild Cartesian axes). It uses the
// exact eigenvalues for Schwarzschild in the body's own frame,
//     (M / r^3) [ (2 + 3K/r^2), -(1 + 3K/r^2), -1 ]
// along the radial, in-plane and normal directions (K = Carter's constant
// Q + (L - aE)^2): stronger than the Newtonian tide for fast, close passes.
Mat3 tidalTensor(double a, const kerr::Constants& c, const kerr::State& s);

// A real star is centrally condensed, so it holds together better than a uniform
// ball of the same mean density: hydrodynamic simulations of gamma = 5/3 stars
// (Guillochon & Ramirez-Ruiz 2013) find full disruption for beta = r_t / r_p
// above about 0.9, where the uniform affine star already breaks up near 0.36.
// Its dynamical frequency is scaled up by this factor (in w^2) to match.
constexpr double kCondensation = 16.0;

// w^2 for a star whose tidal radius is r_t.
inline double frequencySquared(double tidalRadius) {
    return kCondensation / (tidalRadius * tidalRadius * tidalRadius);
}

struct Star {
    Mat3 F = identity();
    Mat3 Fdot{};  // dF/dtau
};

// d^2F/dtau^2 for tidal tensor C and dynamical frequency squared w2.
Mat3 acceleration(const Mat3& F, const Mat3& C, double w2, double gamma = 5.0 / 3.0);

// Advances the star by dtau of its proper time (velocity Verlet, substepped).
void advance(Star& star, const Mat3& C, double w2, double dtau);

// Longest semi-axis in units of the undisturbed radius.
double stretch(const Mat3& F);

// (F F^T)^(-1/2): maps the ellipsoid back onto a sphere of the original radius.
Mat3 inverseShape(const Mat3& F);

}  // namespace tidal
