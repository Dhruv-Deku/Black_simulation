// Shared by the neutron star merger's compute shaders (src/merger_sim.*): the particles,
// the smoothing kernel and the neighbour grid.
//
// Units: G = c = 1 and the Sun's mass is 1 (lengths 1.4766 km, times 4.9255 us); see src/merger.h.

#ifndef MERGER_COMMON
#define MERGER_COMMON

const float PI = 3.14159265358979;

// One parcel of neutron star matter (std430; MergerSim::GpuParticle).
struct Particle {
    vec4 pos;     // xyz, mass
    vec4 vel;     // xyz, smoothing length h
    vec4 acc;     // xyz, d u_th / dt
    vec4 thermo;  // thermal energy per unit mass u_th, density, pressure, star (0 or 1; -1: gone)
    vec4 extra;   // where it started in its star (in stellar radii: follows the matter), potential
};

layout(std430, binding = 0) buffer Particles {
    Particle particles[];
};

// The neighbour grid: a hash table of cells, each the head of a linked list of particles.
layout(std430, binding = 1) buffer GridHead {
    int gridHead[];
};
layout(std430, binding = 2) buffer GridNext {
    int gridNext[];
};

// Per-particle rates the density pass leaves for the force pass: velocity divergence,
// curl magnitude, sound speed.
layout(std430, binding = 4) buffer Aux {
    vec4 aux[];
};

// What the reduction measures each step, and the clock (MergerSim::GpuState).
layout(std430, binding = 3) buffer State {
    double sTime;        // simulated time
    double sPad;
    vec4 sStep;          // dt this step, CFL limit, steps taken, merged (1 once the cores touch)
    vec4 sStar[4];       // star 0: centre of mass, mass; velocity, -; star 1: the same
    vec4 sReaction;      // radiation-reaction acceleration of the relative orbit, on (1) / off (0)
    vec4 sIddot0;        // d^2 I / dt^2: xx, yy, zz, xy
    vec4 sIddot1;        // xz, yz, max density, max temperature (log10 K)
    vec4 sExtent;        // rms radius of the bound matter, outermost dense parcel, unbound mass, kinetic energy
    vec4 sHole;          // black hole: position xyz, mass (0: none)
    vec4 sHoleV;         // its velocity, accreted mass since the last read
    ivec4 sCounters;     // history samples written, particles swallowed, -, -
};

// One sample of the wave history per step (MergerSim::GpuSample).
struct WaveSample {
    double t;
    double pad;
    vec4 a;  // Iddot xx, yy, zz, xy
    vec4 b;  // Iddot xz, yz, separation of the stars' centres, -
};
layout(std430, binding = 5) buffer History {
    WaveSample history[];
};

uniform int uCount;          // particles
uniform float uCellSize;     // neighbour grid cell: about twice the cores' smoothing length
uniform int uHashMask;       // hash table size - 1

// The M4 cubic spline (Monaghan & Lattanzio 1985), support 2h.
float kernelW(float r, float h) {
    float q = r / h;
    float s = 1.0 / (PI * h * h * h);
    if (q < 1.0) return s * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    if (q < 2.0) {
        float t = 2.0 - q;
        return s * 0.25 * t * t * t;
    }
    return 0.0;
}

// dW/dr.
float kernelDW(float r, float h) {
    float q = r / h;
    float s = 1.0 / (PI * h * h * h * h);
    if (q < 1.0) return s * (-3.0 * q + 2.25 * q * q);
    if (q < 2.0) {
        float t = 2.0 - q;
        return -s * 0.75 * t * t;
    }
    return 0.0;
}

ivec3 cellOf(vec3 p) {
    return ivec3(floor(p / uCellSize));
}

int hashCell(ivec3 c) {
    return ((c.x * 73856093) ^ (c.y * 19349663) ^ (c.z * 83492791)) & uHashMask;
}

#endif
