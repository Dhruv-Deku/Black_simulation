#pragma once

#include "central.h"
#include "debris.h"
#include "kerr.h"
#include "penrose.h"
#include "tde.h"
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
    Moon,           // an airless rocky or icy body (moons, Mercury, asteroids, comets), or a satellite
};
const char* kindName(BodyKind kind);

// How a Moon-kind body's surface looks (the GPU's kPalette* in shaders/surfaces.glsl).
enum class Palette { Moon, Io, Europa, Ganymede, Callisto, Mercury, Venus, Comet, Satellite, Rocky, Scorched };

struct Body;
// A representative color for lists, the orbit map and path lines.
glm::vec3 displayColor(const Body& body);

// The radius to draw a body of this radius at this distance with: at least
// kPointRadiusPixels pixels (pixelAngle radians each), so that a body too small to
// resolve still shows as a point; dim is what its brightness is multiplied by so that
// it sends the same total light, (radius / drawn)^2.
constexpr double kPointRadiusPixels = 0.75;
double pointSourceRadius(double radius, double distance, double pixelAngle, double& dim);
// Such points are drawn no fainter than this: a planet seen across its orbit sends a
// pixel's worth far less light than its star does, yet the eye finds it against the night
// sky (Venus outshines every star). The ray tracer's sky is not on a physical scale either.
constexpr double kPointLeastBrightness = 0.25;

// A body following a timelike Kerr geodesic, nudged by the pull of the other bodies.
struct Body {
    std::string name;
    BodyKind kind = BodyKind::Star;
    kerr::Constants constants;
    kerr::State state;
    double mass = 1e-6;           // in units of the hole's mass (a Sun-like star around a million-Sun hole)
    float radius = 0.5f;          // visual radius, in M
    float temperature = 5800.0f;  // surface temperature, K (stars and pulsars)
    // TidalStar: the temperature its gas has inside, on the debris's scale (0: the
    // surface temperature). It sets the debris's pressure once the star is torn apart.
    float interiorTemperature = 0.0f;
    float brightness = 1.0f;
    float seed = 0.0f;            // varies surface patterns
    Palette palette = Palette::Moon;  // Moon kind: which surface
    float spinRate = 0.0f;        // rotation, radians per M of coordinate time
    float tidalRadius = 0.0f;     // TidalStar: where the hole's tide matches the star's gravity
    tidal::Star shape;            // TidalStar: its deformation (the affine model)
    float maxStretch = 1.0f;      // longest axis reached so far, in radii
    bool inside = false;          // has crossed the event horizon
    // The Penrose process (src/penrose.h): the body splits at its closest approach if that
    // is inside this radius (0: never); the fragment with negative energy remembers which
    // split made it (-1: none).
    float splitRadius = 0.0f;
    int penroseSplit = -1;
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
    // 0.5 M a sample (4096 M of history), times the scene scale: big scenes need longer
    // memories. Light from a white dwarf's giant planet, 350 white-dwarf radii out, takes
    // ~2,100 scene units to reach a camera beside the star.
    static constexpr int kColumns = 8192;
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
    TidalDisruptionFlare,
    GrazingStar,
    BinaryStar,
    PlanetWithMoon,
    GasCloud,
    DebrisRing,
    MiniBlackHole,
    WhiteDwarfCompanion,
    // What really orbits a star or a planet, at true sizes and distances.
    Moon,
    SpaceStation,
    GeostationarySatellite,
    AsteroidImpact,
    GalileanMoons,
    CometImpact,
    InnerPlanets,
    HotJupiter,
    SungrazingComet,
    ProximaPlanets,
    GiantInnerPlanets,
    WhiteDwarfPlanet,
    WhiteDwarfAsteroid,
    CataclysmicCompanion,
    InfallingClump,
    Penrose,  // a body splits in the ergosphere: one piece falls in with negative energy
};

struct PresetInfo {
    Preset preset;
    const char* label;
    const char* description;
    bool object;     // listed with the objects rather than the orbit demonstrations
    unsigned hosts;  // the centers it is offered around (central::hostBit masks)
};
const std::vector<PresetInfo>& presetList();
// Whether a preset makes sense around this center.
bool presetAvailable(const PresetInfo& preset, central::Host host);

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
    // What sits at the center (stars, planets, white dwarfs and neutron stars get real
    // companions at real sizes): its kind, physical mass and radius.
    central::Host host = central::Host::BlackHole;
    double centralMassSolar = 1.0;
    double centralRadiusKm = 0.0;
    double centralTemperature = 0.0;  // K
    double centralSpinPeriod = 0.0;   // seconds (0: not turning)
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
    // Every star torn apart so far: when its debris comes back, and what that powers
    // (the flare and the outflow; see tde.h). Newest last.
    std::vector<tde::Event> tidalEvents;
    // Every Penrose split so far, newest last.
    std::vector<penrose::Split> penroseSplits;

    void clear();
    // Sets the scene scale. If it changed a lot, the history's time step follows it and
    // every body's past is laid down again.
    void setScale(double newScale, double spin);
    // How much simulated time a frame covers: when time runs fast the light-delay history
    // samples more coarsely (see retuneHistory).
    void setFrameTime(double dt, double spin);
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
    // and times go out in render units (scale M; see Renderer). A body smaller than a
    // pixel (pixelAngle radians, seen from cameraPosition, Kerr-Schild Cartesian in M) is
    // drawn a pixel wide and dimmed to keep the light it sends: a point, like a star.
    void fillGpu(double spin, bool lightDelay, double cameraR, const glm::dvec3& cameraPosition, double pixelAngle,
                 GpuBodies& out) const;
    // The farthest any body is from the center, in M (0 without bodies).
    double reach() const;

    // Where bodies stop: the horizon, or the surface.
    double boundary(double spin) const;

private:
    friend struct scene::Access;  // a saved scene carries the state below (src/scene_file.cpp)
    bool push(Body body, double spin);
    // The ZAMO-frame velocity at apoapsis of the orbit between rp and ra (in M), tilted
    // out of the equator by tilt. False if there is none.
    bool apoapsisVelocity(double spin, double rp, double ra, double tilt, double v[3]);
    // A body on a real orbit around a star or planet: semi-major axis aKm, eccentricity
    // e, inclination to the central body's equator inclDeg, released at apoapsis;
    // radius in km. Tidally locked (one face toward the center) unless it is given a
    // rotation period (seconds; negative turns backwards).
    bool addRealOrbit(double spin, const std::string& name, BodyKind kind, double radiusKm, double aKm, double e,
                      double inclDeg, double phase, double massRatio, double rotationSeconds = 0.0);
    // Around anything but a black hole the orbit demonstrations carry a small body of
    // a kind that could really be there (plasma near a neutron star, a satellite round a
    // planet, a comet near a star), at its real size.
    void realTestMass(Body& body) const;
    // The presets of real companions (Moon, planets, satellites...).
    void addRealPreset(double spin, Preset preset);
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
    // The Penrose process: the body at `index` splits where it is (src/penrose.h).
    void split(size_t index, double spin);
    // Removes a body that has physically gone, keeping its light in flight.
    void depart(size_t index);
    std::string nextName();

    void retuneHistory(double spin);

    std::mt19937 m_rng{12345};
    int m_counter = 0;
    double m_frameTime = 0.0;  // simulated time per frame, smoothed
};
