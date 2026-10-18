#pragma once

// Where a geodesic is going: exact classification from its constants of motion,
// numerical prediction of its path, and the trace of the light reaching a pixel.

#include "kerr.h"

#include <array>
#include <vector>

namespace trajectory {

using Point = std::array<double, 3>;

// --- Massive particles -----------------------------------------------------------

enum class Fate {
    Bound,     // oscillates between two turning points forever
    Circular,  // sits on a circular orbit
    Plunges,      // reaches the horizon
    Escapes,      // reaches infinity
    HitsSurface,  // reaches the surface of a star or planet
};
const char* fateName(Fate fate);

// The least energy at which a particle with angular momentum L and Carter constant Q
// can be at radius r: R(r) >= 0 exactly when E >= V(r). Plotting V(r) against the
// particle's E shows where it can go: wherever E sits above the curve.
double effectivePotential(double a, double L, double Q, double r);

struct Analysis {
    Fate fate = Fate::Bound;
    bool hasInner = false, hasOuter = false;  // turning points on either side of r
    double rInner = 0.0, rOuter = 0.0;
    bool stableCircular = false;              // for Fate::Circular
    double inclinationDeg = 0.0;              // atan(sqrt(Q) / |L|)
};

// Exact for a test particle: finds the zeros of R(r) nearest the current radius.
// No turning point inward means it falls in; none outward means it escapes. With a
// central body of radius `surface` (0 for a black hole), turning points below the
// surface do not count: the orbit ends there.
Analysis classify(double a, const kerr::Constants& c, const kerr::State& s, double surface = 0.0);

struct Path {
    std::vector<Point> points;  // Cartesian, spaced at least minSpacing apart
    Fate fate = Fate::Bound;    // what happened within the time simulated
    double duration = 0.0;      // coordinate time covered
    double rMin = 0.0, rMax = 0.0;
    int periapses = 0;          // periapsis passages seen
    double precessionDeg = 0.0; // periapsis advance per radial period (if periapses >= 2)
};

// Integrates forward up to maxTime of coordinate time, or until capture (the horizon,
// or the surface of radius `surface` if nonzero) or escape beyond escapeRadius.
Path predict(double a, kerr::Constants c, kerr::State s, double maxTime, double minSpacing, size_t maxPoints = 4000,
             double surface = 0.0, double escapeRadius = 2000.0);

// --- Light rays ------------------------------------------------------------------

struct Sphere {
    Point center;
    double radius;
};

// A small black hole orbiting the big one. No exact spacetime holds two black holes,
// so it lenses in the thin-lens approximation: a ray passing at impact parameter b is
// bent by 4m/b toward it, and is captured if b < 3 sqrt(3) m.
struct Lens {
    Point center;
    double mass;
};

struct RayScene {
    bool disk = true;
    double diskInner = 2.0, diskOuter = 20.0;
    bool diskPrograde = true;
    double escapeRadius = 250.0;
    std::vector<Sphere> spheres;  // orbiting bodies, at their current positions
    std::vector<Lens> lenses;     // orbiting small black holes
    double surface = 0.0;         // radius of the central body's surface; 0 for a black hole
    double surfaceOmega = 0.0;    // its angular velocity d phi / dt
};

enum class RayEnd { Sky, Horizon, Disk, Body, Unfinished, Surface, Lens };
const char* rayEndName(RayEnd end);

// Thin-lens kick by a small black hole of mass m at `center` (Kerr-Schild Cartesian)
// on a photon at its point of closest approach, (c, s), passing at distance b: its
// direction, measured by the local ZAMO, turns by 4m/b toward the lens at fixed local
// frequency, and E, L, Q and the velocities are recomputed from the new direction.
// With `backward` (tracing from the camera into the past) the photon holds its later
// direction and gets its earlier one, turned the other way. Returns false (and leaves
// the photon alone) if b < 3 sqrt(3) m: the photon is captured.
bool lensKick(double a, const Point& center, double mass, bool backward, kerr::Constants& c, kerr::State& s);

struct RayReport {
    RayEnd end = RayEnd::Unfinished;
    double energy = 1.0;      // conserved energy E of the photon, for unit energy at the camera
    double b = 0.0, q = 0.0;  // impact parameter L/E and Carter constant Q/E^2
    double gCamera = 1.0;     // frequency the camera measures for p_t = -1 (1 / energy)
    double g = 1.0;           // observed / emitted frequency at the source (disk, body or surface)
    double emissionCos = 0.0; // surface: cosine of the emission angle from the normal, in the surface's frame
    double surfaceTheta = 0.0, surfacePhi = 0.0;  // surface: where it left (Kerr-Schild angles)
    double rEnd = 0.0;        // where it was emitted (disk or body), or where it was lost
    int body = -1;
    double delay = 0.0;       // coordinate time the light spent travelling
    double rMin = 0.0;        // closest approach to the hole
    double turns = 0.0;       // how many times it went around the hole
    int lensings = 0;         // deflections by small black holes
    Point skyDirection{};     // for RayEnd::Sky
    std::vector<Point> path;  // from the camera backwards along the ray
};

// The light arriving at a camera at (r, theta, phi) with the given frame (see
// kerr::observerFrame) from direction `look`, traced back to its source. The camera
// may be inside the horizon.
RayReport traceRay(double a, double r, double theta, double phi, const kerr::Frame& frame, const double look[3],
                   const RayScene& scene);

}  // namespace trajectory
