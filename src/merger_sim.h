#pragma once

#include "gl_util.h"
#include "merger.h"

#include <glm/glm.hpp>

#include <array>
#include <filesystem>
#include <string>
#include <vector>

// Two neutron stars as smoothed-particle hydrodynamics on the GPU (shaders/merger_*.comp):
// self-gravity, pressure (a cold n = 1 polytrope plus heat), shocks, the radiation reaction
// that drives the inspiral, and the black hole the remnant collapses into. src/merger.h
// has the physics and its references.
//
// The run goes through these phases:
//   Inspiral   the stars, built and settled into equilibrium, are carried whole along
//              Peters' inspiral: while they are far apart, tides do nothing worth the cost;
//   Hydro      from a few stellar radii apart every parcel is followed: tides stretch the
//              stars, they merge, and the remnant spins and rings;
//   Collapsed  at the chosen moment the remnant collapses into a black hole and what has
//              too much angular momentum stays out as a torus;
//   Coasting   the hydrodynamics stops; the ejecta coast outward (drawn where they would be).
//
// Units: G = c = 1 and the Sun's mass is 1 (src/merger.h).
class MergerSim {
public:
    struct Params {
        double m1 = 1.35, m2 = 1.35;   // solar masses (GW170817-like)
        double radiusKm = 12.0;        // both stars (an n = 1 polytrope's radius does not depend on its mass)
        int perStar = 12000;           // particles in the heavier star (the lighter has as many, lighter ones)
        double startKm = 100.0;        // separation at the start
        double hydroKm = 40.0;         // separation where the hydrodynamics takes over
        double collapseMs = 15.0;      // the remnant collapses this long after the merger
        double coastMs = 8.0;          // the hydrodynamics runs this long after the collapse
        double initialTemperature = 1e8;  // K, the stars' interiors (cold, for neutron stars)
    };
    enum class Phase { Inspiral, Hydro, Collapsed, Coasting };

    bool init(const std::filesystem::path& shaderDir, std::string& error);
    void shutdown();
    void hotReload();
    std::string errors() const;

    // Builds and settles the stars (a second or so of GPU time).
    void reset(const Params& params);
    bool ready() const { return m_count > 0; }
    const Params& params() const { return m_params; }

    // Runs until `target` (units) since the start, or until the step budget is spent. The
    // time actually reached is time().
    void advance(double target, int maxSteps);

    // Fills the drawing grids: a cube of half-size `half` around `center`, and a map from
    // above of half-size mapHalf around mapCenter.
    void deposit(const glm::vec3& center, float half, const glm::vec2& mapCenter, float mapHalf);
    GLuint volumeA() const { return m_volumeA; }
    GLuint volumeB() const { return m_volumeB; }
    GLuint map() const { return m_map; }
    static constexpr int kVolumeDims = 128;
    static constexpr int kMapDims = 1024;
    glm::vec3 boxMin() const { return m_boxMin; }
    float boxSize() const { return m_boxSize; }
    glm::vec2 mapMin() const { return m_mapMin; }
    float mapSize() const { return m_mapSize; }

    // --- What is happening ---------------------------------------------------------------
    Phase phase() const { return m_phase; }
    double time() const { return m_time; }                  // since the start
    // Since the merger (negative before: predicted from Peters' law until it happens).
    double sinceMerger() const;
    bool merged() const { return m_mergeTime >= 0.0; }
    double mergeTime() const { return m_mergeTime; }
    double collapseTime() const { return m_mergeTime >= 0.0 ? m_mergeTime + m_params.collapseMs * 1e-3 / merger::kSecond : 1e30; }
    double separation() const { return m_separation; }      // of the stars' centres
    double radius() const { return m_radius; }              // the stars' radius (units)
    glm::dvec3 starCenter(int k) const { return m_center[k]; }
    double starMass(int k) const { return k == 0 ? m_params.m1 : m_params.m2; }
    double totalMass() const { return m_params.m1 + m_params.m2; }
    double maxDensity() const { return m_maxDensity; }      // g/cm^3
    double maxTemperature() const { return m_maxTemperature; }  // K
    double extent() const { return m_extent; }               // outermost dense matter, from the centre
    double rmsRadius() const { return m_rms; }               // of the bound matter
    double unboundMass() const { return m_unbound; }         // solar masses
    double aliveMass() const { return m_aliveMass; }
    const merger::BlackHole& hole() const { return m_hole; }
    bool hasHole() const { return m_hole.mass > 0.0; }
    double holeMass() const;                                 // grows as it swallows gas
    double gwFrequency() const { return m_gwFrequency; }    // cycles per unit time, from the wave's phase
    double lastStep() const { return m_lastDt; }
    int stepsTaken() const { return m_steps; }
    float gpuMilliseconds() const { return m_gpuMs; }

    // The wave history: time and d^2 I / dt^2 (xx, yy, zz, xy, xz, yz), one entry per step.
    struct WavePoint {
        double t;
        std::array<double, 6> iddot;
    };
    const std::vector<WavePoint>& waves() const { return m_waves; }
    // The strain (+, x) seen from direction n at distance D a time t (units since the
    // start), with the hole's ringdown added after the collapse.
    merger::Strain strainAt(double t, const glm::dvec3& n, double distance) const;
    // d^2 I / dt^2 (xx, yy, zz, xy, xz, yz) at time t, ringdown included.
    std::array<double, 6> iddotAt(double t) const;

    // Unbound matter, for the kilonova: mass-weighted directions of its velocities (a
    // histogram over the sky, 32 x 16 bins in azimuth and polar angle) and its mean speed.
    const std::vector<float>& ejectaSky() const { return m_ejectaSky; }
    double ejectaSpeed() const { return m_ejectaSpeed; }

private:
    struct GpuParticle {
        glm::vec4 pos, vel, acc, thermo, extra;
    };
    struct GpuState {
        double time;
        double pad;
        glm::vec4 step;
        glm::vec4 star[4];
        glm::vec4 reaction;
        glm::vec4 iddot0, iddot1;
        glm::vec4 extent;
        glm::vec4 hole, holeV;
        glm::ivec4 counters;
    };
    struct GpuSample {
        double t;
        double pad;
        glm::vec4 a, b;
    };
    static_assert(sizeof(GpuState) == 16 + 16 + 64 + 16 * 6 + 16, "std430 layout");
    static_assert(sizeof(GpuSample) == 48, "std430 layout");

    void bindAll();
    void setCommon(ShaderProgram& p);
    void dispatchStep(double dtMax, bool record);
    void computeForces();
    void readState();
    void carry(double t);  // the stars' centres on Peters' inspiral at time t
    void startHydro();
    void doCollapse();
    void measureEjecta();
    std::vector<GpuParticle> readParticles() const;
    void sortParticles();
    int m_sortedAt = -1000;

    Params m_params;
    ShaderProgram m_grid, m_density, m_force, m_integrate, m_reduce, m_deposit, m_resolve;
    GLuint m_particles = 0, m_head = 0, m_next = 0, m_state = 0, m_aux = 0, m_history = 0, m_cells = 0, m_counts = 0;
    static constexpr int kCountDims = 48;  // the coarse grid that measures how sparse the gas is
    GLuint m_volumeA = 0, m_volumeB = 0, m_map = 0;
    GLuint m_query = 0;
    int m_count = 0;
    int m_capacity = 0;
    static constexpr int kHashSize = 1 << 19;
    static constexpr int kHistorySize = 1 << 16;
    int m_historyRead = 0;

    double m_radius = 8.0;         // stars' radius, units
    double m_K = 40.0;
    double m_hCore = 0.5;
    double m_cell = 1.0;
    double m_massRef = 1e-4;

    Phase m_phase = Phase::Inspiral;
    double m_time = 0.0;
    double m_mergeTime = -1.0;
    double m_coastStart = -1.0;
    double m_lastDt = 0.5;
    int m_steps = 0;
    // Analytic inspiral (Inspiral phase): separation and orbital phase every m_orbitDt from
    // the start to m_hydroStart, when the hydrodynamics takes over.
    std::vector<glm::dvec2> m_orbit;
    double m_orbitDt = 1.0;
    double m_hydroStart = 0.0;
    double m_mergeSeparation = 10.0;
    std::array<glm::dvec3, 2> m_center{}, m_velocity{};
    double m_separation = 0.0;
    double m_maxDensity = 0.0, m_maxTemperature = 0.0, m_extent = 0.0, m_rms = 0.0, m_unbound = 0.0;
    double m_aliveMass = 0.0;
    merger::BlackHole m_hole;
    double m_holeTime = -1.0;
    merger::Ringdown m_ring;
    double m_ringAmplitude[2] = {0.0, 0.0};  // the remnant's quadrupole at collapse, to ring down from
    std::vector<WavePoint> m_waves;
    double m_gwFrequency = 0.0;
    std::vector<float> m_ejectaSky;
    double m_ejectaSpeed = 0.2;
    glm::vec3 m_boxMin{0.0f};
    float m_boxSize = 1.0f;
    glm::vec2 m_mapMin{0.0f};
    float m_mapSize = 1.0f;
    float m_gpuMs = 0.0f;
};
