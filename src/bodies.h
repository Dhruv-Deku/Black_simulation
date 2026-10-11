#pragma once

#include "debris.h"
#include "kerr.h"
#include "tidal.h"

#include <glm/glm.hpp>

#include <deque>
#include <random>
#include <string>
#include <vector>

enum class BodyKind {
    Star,
    RockyPlanet,
    GasGiant,
    Pulsar,
    TidalStar,      // a star whose shape responds to tides, and is torn apart by them
    WhiteDwarf,
    NeutronStar,
    MiniBlackHole,  // lenses light passing it (thin-lens approximation) and swallows what comes too close
};
const char* kindName(BodyKind kind);

struct Body;
// A representative color for lists, the orbit map and path lines.
glm::vec3 displayColor(const Body& body);

// A body following a timelike Kerr geodesic, nudged by the pull of the other bodies.
struct Body {
    std::string name;
    BodyKind kind = BodyKind::Star;
    kerr::Constants constants;
    kerr::State state;
    double mass = 1e-6;           // in units of the hole's mass (a Sun-like star around a million-Sun hole)
    float radius = 0.5f;          // visual radius, in M
    float temperature = 5800.0f;  // surface temperature, K (stars and pulsars)
    float brightness = 1.0f;
    float seed = 0.0f;            // varies surface patterns
    float spinRate = 0.0f;        // rotation, radians per M of coordinate time
    float tidalRadius = 0.0f;     // TidalStar: where the hole's tide matches the star's gravity
    tidal::Star shape;            // TidalStar: its deformation (the affine model)
    float maxStretch = 1.0f;      // longest axis reached so far, in radii
    bool inside = false;          // has crossed the event horizon
    int slot = -1;                // row in the light-delay history
    std::deque<glm::vec3> trail;  // recent Cartesian positions, for the orbit map
};

// Past positions, 4-velocities and shapes of every body, sampled every `step` of
// coordinate time, so the ray tracer can draw each body where (and as) it was when
// the light now arriving at the camera left it. A ring buffer of kColumns samples
// per body, mirrored in a GPU texture with kRows rows per slot:
//   0: x, y, z, r    1: u^t, u^r, u^theta, u^phi
//   2: S_xx, S_yy, S_zz, stretch    3: S_xy, S_xz, S_yz, 0
// where S = (F F^T)^(-1/2) maps the (tidally deformed) body back onto a sphere.
struct History {
    static constexpr int kSlots = 32;
    static constexpr int kRows = 4;
    static constexpr int kColumns = 4096;
    // 0.5 M (2048 M of history), times the scene scale: big scenes need longer memories.
    double step = 0.5;

    struct Sample {
        glm::vec4 position, u, shapeA, shapeB;
    };

    std::vector<glm::vec4> image = std::vector<glm::vec4>(static_cast<size_t>(kColumns) * kRows * kSlots);
    long long newest = 0;       // index of the newest sample; its time is newest * step
    int length[kSlots] = {};    // valid samples per slot
    std::vector<int> dirtyRows; // slots whose whole row must be re-uploaded
    int generation = 0;         // bumped when the step changes: everything must be re-uploaded

    void record(int slot, long long k, const Sample& sample);
};

// std140 mirror of the Bodies uniform block in shaders/trace.frag.
struct GpuBody {
    glm::vec4 positionRadius;  // current Cartesian center, radius
    glm::vec4 fourVelocity;    // current u^t, u^r, u^theta, u^phi
    glm::vec4 params;          // temperature, brightness, kind, seed
    glm::vec4 motion;          // rotation angle now, rotation rate, history slot, history length
    glm::vec4 range;           // r range over the light-delay window (min, max), time since lost (-1 if not),
                               // bounding radius
    glm::vec4 shapeA;          // current shape, as in History rows 2 and 3
    glm::vec4 shapeB;
};
struct GpuBodies {
    static constexpr int kMax = History::kSlots;
    GpuBody bodies[kMax];
    glm::ivec4 info;  // count, history head column, light delay on, history columns
    glm::vec4 shell;  // r range of all bodies (min, max), fraction of a step since the newest sample, step
};
static_assert(sizeof(GpuBody) == 112, "std140 layout");
static_assert(sizeof(GpuBodies) == 112 * 32 + 32, "std140 layout");

enum class Preset {
    IscoOrbit,
    StableOrbit,
    RetrogradeOrbit,
    Precessing,
    ZoomWhirl,
    PolarOrbit,
    Plunge,
    Flyby,
    Cluster,
    RockyPlanet,
    GasGiant,
    Pulsar,
    TidalDisruption,
    GrazingStar,
    BinaryStar,
    PlanetWithMoon,
    GasCloud,
    DebrisRing,
    MiniBlackHole,
    WhiteDwarfCompanion,
};

struct PresetInfo {
    Preset preset;
    const char* label;
    const char* description;
    bool object;  // listed with the objects rather than the orbit demonstrations
};
const std::vector<PresetInfo>& presetList();

// A body that struck the surface of the central star or planet.
struct Impact {
    double theta, phi;  // where (Kerr-Schild angles)
    double time;        // when (coordinate time)
    double energy;      // released per unit rest mass, (gamma - 1) in the surface's frame
    double mass;        // in units of M
};

class BodySystem {
public:
    double time = 0.0;  // coordinate time t, in M
    std::vector<Body> bodies;
    // The central object, as the bodies need it (kept up to date by the app):
    double surface = 0.0;       // surface radius in M; 0 for a black hole (bodies stop at the horizon)
    double surfaceOmega = 0.0;  // its angular velocity
    double scale = 1.0;         // scene scale (see Settings): sizes in presets are multiples of it
    double massSolar = 0.0;     // physical mass of the hole when physical units are on (true-scale presets)
    std::vector<Impact> impacts;  // recent impacts on the surface, oldest first (the app trims them)
    double accreted = 0.0;        // mass that has hit the surface so far, in M
    int impactCount = 0;          // bodies that have hit the surface so far
    // Bodies that are gone (torn apart, swallowed, escaped) but whose last light may
    // still be on its way: drawn from their history until it has all arrived.
    struct Departed {
        Body body;
        double lostTime;
    };
    std::vector<Departed> departed;
    std::deque<std::string> events;              // recent happenings, newest first
    History history;
    // Particles waiting to be handed to the GPU. Gas from a preset (a cloud, a ring) has
    // been there all along; debris from a disruption only exists from now on.
    struct DebrisBatch {
        std::vector<ParticleInit> particles;
        bool existedBefore = false;
    };
    std::vector<DebrisBatch> pendingDebris;

    void clear();
    // Sets the scene scale. If it changed a lot, the history's time step follows it and
    // every body's past is laid down again.
    void setScale(double newScale, double spin);
    // Advances everything by dt. With gravity on, bodies also pull on each other.
    void update(double spin, double dt, bool gravity);
    // Keeps every body's position and locally measured velocity when the spin changes.
    void respin(double oldSpin, double newSpin);

    // Each returns false (and logs why) if the body could not be created.
    bool addFromZamo(double spin, double r, double theta, double phi, const double v[3],
                     BodyKind kind = BodyKind::Star, const std::string& name = {});
    bool addCircular(double spin, double r, bool prograde, double phase, BodyKind kind = BodyKind::Star,
                     const std::string& name = {});
    bool addEquatorial(double spin, double rp, double ra, bool prograde, double phase,
                       BodyKind kind = BodyKind::Star, const std::string& name = {});
    void addPreset(double spin, Preset preset);
    // Removes a body at once (the user deleted it).
    void remove(size_t index);
    void logEvent(std::string text);

    // cameraR sets how far back in time the light-delay culling range must reach. Lengths
    // and times go out in render units (scale M; see Renderer).
    void fillGpu(double spin, bool lightDelay, double cameraR, GpuBodies& out) const;

    // Where bodies stop: the horizon, or the surface.
    double boundary(double spin) const;

private:
    bool push(Body body, double spin);
    // Small black holes swallow bodies that come within their horizon (plus the body's size).
    void swallowByMiniHoles(double spin);
    void backfill(Body& body, double spin);
    void recordTrail(Body& body, double spin);
    // Moves every body along its geodesic to coordinate time `target`.
    void drift(double spin, double target);
    // Adds a body of the given mass orbiting another one (a moon, or a binary companion).
    bool addCompanion(double spin, size_t primary, double mass, double separation, double tiltDeg, BodyKind kind,
                      const std::string& name);
    void disrupt(const Body& body, double spin);
    // Removes a body that has physically gone, keeping its light in flight.
    void depart(size_t index);
    std::string nextName();

    std::mt19937 m_rng{12345};
    int m_counter = 0;
};
