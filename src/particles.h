#pragma once

#include "debris.h"
#include "gl_util.h"

#include <glm/glm.hpp>

#include <filesystem>
#include <string>
#include <vector>

// How the debris gas behaves (see shaders/particles_update.comp).
struct GasParams {
    bool physics = true;             // pressure, shocks, heating and cooling
    float collisionRate = 0.3f;      // per M, per unit density: how readily crossing streams shock
    float coolingTime = 60.0f;       // M, at 10,000 K (hotter gas cools faster, ~T^-3)
    float floorTemperature = 4500.0f;
    std::vector<glm::vec4> bodies;   // Cartesian position and mass of each body pulling on the gas (M)

    // The central object: scene scale (the gas lives in render units of scale M), and
    // for a star or planet its surface (M), where parcels stop and leave hot spots.
    double scale = 1.0;
    double surface = 0.0;
    double surfaceShift = 0.0;   // phi - phi_BL at the surface
    double surfacePhase = 0.0;   // rotation angle of the surface now
    double impactFadeTime = 20.0;  // e-folding time of the hot spots, in M

    // Light delay: draw each parcel where it was when the light now reaching the
    // camera (at this Cartesian position and radius, outside the horizon) left it.
    bool lightDelay = false;
    glm::vec3 cameraPosition{0.0f};
    float cameraR = 0.0f;
    double time = 0.0;  // coordinate time after this update (M)
};

// Debris on the GPU. Each particle is a parcel of gas following its own Kerr
// geodesic (compute shader), nudged by the gas around it and by the bodies. Each
// frame the parcels are deposited into a (log r, theta, phi) grid that the physics
// reads back the next frame (full precision, ZAMO-frame velocities), and into a
// second grid that the ray tracer samples as glowing, absorbing gas (with light
// delay, their light-delayed copies), so the debris is lensed and Doppler shifted
// like everything else. On the GPU everything is in render units of scale() M (see
// particle_common.glsl); the interface takes M.
class ParticleSystem {
public:
    static constexpr int kCapacity = 1 << 19;  // 524,288 particles
    // Each spawned object carries this much "mass", however many particles it uses,
    // so its brightness does not depend on the particle count.
    static constexpr float kObjectMass = 100000.0f;
    static constexpr int kGridR = 160, kGridTheta = 96, kGridPhi = 256;
    static constexpr float kGridRMin = 1.0f, kGridRMax = 100.0f;

    bool init(const std::filesystem::path& shaderDir, std::string& error);
    void shutdown();

    // Adds particles at coordinate time `time`, overwriting the oldest once the capacity
    // is reached. If they existedBefore, their past is seen at once (followed back along
    // their orbits); otherwise they appear once the light from their birth arrives.
    void spawn(const std::vector<ParticleInit>& particles, double time, bool existedBefore);
    void clear();
    // The spin (a, in M) or the scene scale changed: keeps every parcel's place and local
    // velocity. Parcels inside the new horizon or surface (M) are dropped.
    void reframe(double oldSpin, double newSpin, double oldScale, double newScale, double surface);
    // Advances every particle by dt of coordinate time and rebuilds the grids.
    void update(double spin, double dt, const GasParams& gas);

    bool active() const { return m_used > 0; }
    int alive() const { return m_alive; }
    // What the ray tracer draws.
    GLuint densityTexture() const { return m_drawDensity; }
    GLuint velocityTexture() const { return m_drawVelocity; }
    GLuint occupancyTexture() const { return m_occupancy; }
    // Hot spots left by gas striking a surface (0 until any has).
    GLuint impactTexture() const { return m_impactMap; }
    static constexpr int kImpactTheta = 64, kImpactLongitude = 128;
    // Weight of gas that has struck the surface so far (1 = one whole cloud or star).
    double accreted() const { return m_accreted; }
    double scale() const { return m_scale; }
    static constexpr int kOccupancyR = kGridR / 8, kOccupancyTheta = kGridTheta / 8, kOccupancyPhi = kGridPhi / 8;

    // A summary of the particles (reads them back from the GPU: slow, for diagnostics).
    std::string describe() const;

    void hotReload();
    std::string errors() const;

private:
    void allocate();
    void resolve(GLuint density, GLuint velocity, GLenum format, bool physics);

    ShaderProgram m_update, m_deposit, m_resolve, m_respin, m_impactResolve;
    GLuint m_particles = 0;  // SSBO of particles
    GLuint m_grid = 0;       // SSBO of fixed-point voxel sums
    GLuint m_stats = 0;      // SSBO with the alive count
    GLuint m_density = 0;    // 3D textures of the gas as it is now, read back by the physics
    GLuint m_velocity = 0;
    GLuint m_drawDensity = 0;  // ... and as the ray tracer draws it
    GLuint m_drawVelocity = 0;
    GLuint m_occupancy = 0;  // coarse map of where there is any debris to draw
    GLuint m_impacts = 0;    // SSBO: energy left on the surface this frame
    GLuint m_impactMap = 0;  // ... accumulated and fading (R32F, colatitude x longitude)
    double m_scale = 1.0;    // render unit in M
    double m_accreted = 0.0;
    int m_used = 0;  // slots holding particles
    int m_next = 0;  // next slot to write
    int m_alive = 0;
    int m_frame = 0;
    bool m_dirty = false;  // grid needs rebuilding even if no time passed
};
