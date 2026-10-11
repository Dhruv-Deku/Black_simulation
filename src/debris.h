#pragma once

// Initial conditions for debris particles: tidally disrupted stars, gas clouds and
// rings. Each particle is a test mass on its own Kerr geodesic; the GPU integrates
// them (shaders/particles_update.comp).

#include "kerr.h"
#include "tidal.h"

#include <random>
#include <vector>

struct ParticleInit {
    float r, theta, phi;
    float vr, vtheta;  // Mino-time velocities
    float E, L, Q;
    float temperature;  // K
};

namespace debris {

// A cloud of `count` particles spread through a ball of the given radius around
// (r, theta, phi) (uniformly, or Gaussian with that standard deviation). Every
// particle starts with the same velocity *vector* as the center, v being the
// center's ZAMO-frame velocity: the spread in position becomes a spread in orbital
// energy that shears the cloud into a stream. Particles that would start inside the
// central body (radius `surface`; 0 for a black hole: inside the horizon) are left out.
std::vector<ParticleInit> cloud(double a, double r, double theta, double phi, const double v[3], double radius,
                                bool gaussian, int count, float temperature, std::mt19937& rng,
                                double surface = 0.0);

// The debris of a star of the given radius torn apart while deformed into
// x = F X (see tidal.h), with internal velocities dF/dtau X on top of the center's
// ZAMO-frame velocity v. F and Fdot are on Boyer-Lindquist Cartesian axes.
std::vector<ParticleInit> deformedStar(double a, double r, double theta, double phi, const double v[3],
                                       double radius, const tidal::Mat3& F, const tidal::Mat3& Fdot, int count,
                                       float temperature, std::mt19937& rng, double surface = 0.0);

// A thin ring of particles on near-circular orbits between rInner and rOuter,
// tilted by tiltDeg about the x axis. Frame dragging precesses tilted orbits at a
// rate that falls off as r^-3, so the ring slowly twists into a warp.
std::vector<ParticleInit> ring(double a, double rInner, double rOuter, double tiltDeg, int count,
                               float temperature, std::mt19937& rng, double surface = 0.0);

}  // namespace debris
