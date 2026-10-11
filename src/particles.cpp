#include "particles.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {

struct GpuParticle {
    glm::vec4 state;         // r, theta, phi, life (1 alive, 0.5 gone but still seen, 0 dead)
    glm::vec4 velocity;      // vr, vtheta, temperature, weight
    glm::vec4 constants;     // E, L, Q, gas density last frame
    glm::vec4 seen;          // the light-delayed copy: r, theta, phi, and the time it is at
    glm::vec4 seenVelocity;  // vr, vtheta, time created (-1e30: always existed), time lost
};
static_assert(sizeof(GpuParticle) == 80, "std430 layout");

constexpr int kChannels = 6;
constexpr float kDensityScale = 1e-3f;
constexpr float kKelvinPerC2 = 1.3e7f;  // see particles_update.comp

glm::ivec3 gridDims() {
    return {ParticleSystem::kGridR, ParticleSystem::kGridTheta, ParticleSystem::kGridPhi};
}

glm::vec2 gridLogR() {
    return {std::log(ParticleSystem::kGridRMin), std::log(ParticleSystem::kGridRMax)};
}

}  // namespace

bool ParticleSystem::init(const std::filesystem::path& shaderDir, std::string& error) {
    bool ok = m_update.buildCompute(shaderDir, "particles_update.comp");
    ok &= m_deposit.buildCompute(shaderDir, "particles_deposit.comp");
    ok &= m_resolve.buildCompute(shaderDir, "grid_resolve.comp");
    ok &= m_respin.buildCompute(shaderDir, "particles_respin.comp");
    ok &= m_impactResolve.buildCompute(shaderDir, "impact_resolve.comp");
    if (!ok) error = errors();
    return ok;
}

void ParticleSystem::allocate() {
    if (m_particles) return;
    glCreateBuffers(1, &m_particles);
    glNamedBufferStorage(m_particles, sizeof(GpuParticle) * kCapacity, nullptr, GL_DYNAMIC_STORAGE_BIT);
    glClearNamedBufferData(m_particles, GL_R32F, GL_RED, GL_FLOAT, nullptr);

    const size_t voxels = static_cast<size_t>(kGridR) * kGridTheta * kGridPhi;
    glCreateBuffers(1, &m_grid);
    glNamedBufferStorage(m_grid, voxels * kChannels * sizeof(int), nullptr, GL_DYNAMIC_STORAGE_BIT);

    glCreateBuffers(1, &m_stats);
    glNamedBufferStorage(m_stats, 2 * sizeof(int), nullptr, GL_DYNAMIC_STORAGE_BIT);

    glCreateBuffers(1, &m_impacts);
    glNamedBufferStorage(m_impacts, sizeof(int) * kImpactTheta * kImpactLongitude, nullptr, GL_DYNAMIC_STORAGE_BIT);
    glClearNamedBufferData(m_impacts, GL_R32I, GL_RED_INTEGER, GL_INT, nullptr);
    glCreateTextures(GL_TEXTURE_2D, 1, &m_impactMap);
    glTextureStorage2D(m_impactMap, 1, GL_R32F, kImpactTheta, kImpactLongitude);
    glTextureParameteri(m_impactMap, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(m_impactMap, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(m_impactMap, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);  // colatitude
    glTextureParameteri(m_impactMap, GL_TEXTURE_WRAP_T, GL_REPEAT);         // longitude
    glClearTexImage(m_impactMap, 0, GL_RED, GL_FLOAT, nullptr);

    auto makeVolume = [](GLuint& tex, GLenum format) {
        glCreateTextures(GL_TEXTURE_3D, 1, &tex);
        glTextureStorage3D(tex, 1, format, kGridR, kGridTheta, kGridPhi);
        glTextureParameteri(tex, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTextureParameteri(tex, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTextureParameteri(tex, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);  // log r
        glTextureParameteri(tex, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);  // theta
        glTextureParameteri(tex, GL_TEXTURE_WRAP_R, GL_REPEAT);         // phi wraps around
        // The first update reads the gas as it was "last frame" before anything is written.
        glClearTexImage(tex, 0, GL_RGBA, GL_FLOAT, nullptr);
    };
    // The physics needs full precision: relative velocities of a few percent of the
    // flow decide how much the gas heats up.
    makeVolume(m_density, GL_RGBA32F);
    makeVolume(m_velocity, GL_RGBA32F);
    makeVolume(m_drawDensity, GL_RGBA16F);
    makeVolume(m_drawVelocity, GL_RGBA16F);

    // (Building this with glGenerateTextureMipmap instead took over 100 ms a frame.)
    glCreateTextures(GL_TEXTURE_3D, 1, &m_occupancy);
    glTextureStorage3D(m_occupancy, 1, GL_R32UI, kOccupancyR, kOccupancyTheta, kOccupancyPhi);
    glTextureParameteri(m_occupancy, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(m_occupancy, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}

void ParticleSystem::shutdown() {
    GLuint buffers[] = {m_particles, m_grid, m_stats, m_impacts};
    glDeleteBuffers(4, buffers);
    GLuint textures[] = {m_density, m_velocity, m_drawDensity, m_drawVelocity, m_occupancy, m_impactMap};
    glDeleteTextures(6, textures);
    m_particles = m_grid = m_stats = m_impacts = 0;
    m_density = m_velocity = m_drawDensity = m_drawVelocity = m_occupancy = m_impactMap = 0;
    m_used = m_next = m_alive = 0;
}

void ParticleSystem::spawn(const std::vector<ParticleInit>& particles, double time, bool existedBefore) {
    if (particles.empty()) return;
    allocate();
    const float weight = kObjectMass / static_cast<float>(particles.size());
    // Into render units: lengths, L and v_theta over the scale; Q and v_r over its
    // square; times over the scale.
    const double k = 1.0 / m_scale;
    const float now = static_cast<float>(time * k);
    const float created = existedBefore ? -1e30f : now;
    std::vector<GpuParticle> data;
    data.reserve(particles.size());
    for (const ParticleInit& p : particles) {
        const float r = static_cast<float>(p.r * k);
        const float vr = static_cast<float>(p.vr * k * k), vth = static_cast<float>(p.vtheta * k);
        data.push_back({glm::vec4(r, p.theta, p.phi, 1.0f), glm::vec4(vr, vth, p.temperature, weight),
                        glm::vec4(p.E, static_cast<float>(p.L * k), static_cast<float>(p.Q * k * k), 0.0f),
                        glm::vec4(r, p.theta, p.phi, now), glm::vec4(vr, vth, created, 1e30f)});
    }
    // Write into the ring of slots, wrapping around over the oldest particles.
    size_t done = 0;
    while (done < data.size()) {
        size_t n = std::min(data.size() - done, static_cast<size_t>(kCapacity - m_next));
        glNamedBufferSubData(m_particles, sizeof(GpuParticle) * m_next, sizeof(GpuParticle) * n, data.data() + done);
        done += n;
        m_next = (m_next + static_cast<int>(n)) % kCapacity;
        m_used = std::min(kCapacity, m_used + static_cast<int>(n));
    }
    m_alive += static_cast<int>(data.size());
    m_dirty = true;
}

void ParticleSystem::clear() {
    m_used = m_next = m_alive = 0;
}

void ParticleSystem::reframe(double oldSpin, double newSpin, double oldScale, double newScale, double surface) {
    m_scale = newScale;
    if (!active() || !m_respin.valid()) return;
    m_respin.use();
    m_respin.set("uCount", m_used);
    m_respin.set("uRescale", static_cast<float>(oldScale / newScale));
    m_respin.set("uOldSpin", static_cast<float>(oldSpin / newScale));
    m_respin.set("uSpin", static_cast<float>(newSpin / newScale));
    m_respin.set("uMass", static_cast<float>(1.0 / newScale));
    m_respin.set("uSurfaceR", static_cast<float>(surface / newScale));
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_particles);
    glDispatchCompute((m_used + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    m_dirty = true;
}

void ParticleSystem::resolve(GLuint density, GLuint velocity, GLenum format, bool physics) {
    m_resolve.use();
    m_resolve.set("uGridDims", gridDims());
    m_resolve.set("uGridLogR", gridLogR());
    m_resolve.set("uDensityScale", kDensityScale);
    m_resolve.set("uOccupancyDims", glm::ivec3(kOccupancyR, kOccupancyTheta, kOccupancyPhi));
    m_resolve.set("uMarkOccupancy", physics ? 0 : 1);
    m_resolve.set("uPhysics", physics ? 1 : 0);
    glBindImageTexture(0, density, 0, GL_TRUE, 0, GL_WRITE_ONLY, format);
    glBindImageTexture(1, velocity, 0, GL_TRUE, 0, GL_WRITE_ONLY, format);
    glBindImageTexture(2, m_occupancy, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R32UI);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_grid);
    glDispatchCompute((kGridR + 7) / 8, (kGridTheta + 7) / 8, (kGridPhi + 3) / 4);
}

void ParticleSystem::update(double spin, double dt, const GasParams& gas) {
    if (!active() || !m_update.valid() || !m_resolve.valid() || !m_deposit.valid()) return;
    if (dt <= 0.0 && !m_dirty) return;
    m_dirty = false;
    const bool delayed = gas.lightDelay;

    glClearNamedBufferData(m_grid, GL_R32I, GL_RED_INTEGER, GL_INT, nullptr);
    // Keep the accreted count (second int); clear only the alive count.
    const int zeroInt = 0;
    glNamedBufferSubData(m_stats, 0, sizeof(int), &zeroInt);
    const GLuint zero = 0;
    glClearTexImage(m_occupancy, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &zero);

    // Render units: the scene scale is the length unit, and the central mass is 1 / scale.
    m_scale = gas.scale;
    const double k = 1.0 / m_scale;
    const double dynamical = std::sqrt(m_scale);  // dynamical times grow as sqrt(scale) in these units
    m_update.use();
    m_update.set("uCount", m_used);
    m_update.set("uSpin", static_cast<float>(spin * k));
    m_update.set("uMass", static_cast<float>(k));
    m_update.set("uSurfaceR", static_cast<float>(gas.surface * k));
    m_update.set("uKelvinPerC2", static_cast<float>(kKelvinPerC2 * m_scale));
    m_update.set("uDt", static_cast<float>(dt * k));
    m_update.set("uGridDims", gridDims());
    m_update.set("uGridLogR", gridLogR());
    m_update.set("uGasPhysics", gas.physics ? 1 : 0);
    m_update.set("uCollisionRate", static_cast<float>(gas.collisionRate / dynamical));
    m_update.set("uCoolingTime", static_cast<float>(std::max(gas.coolingTime, 1.0f) * dynamical));
    m_update.set("uFloorTemperature", gas.floorTemperature);
    const int bodyCount = static_cast<int>(std::min<size_t>(gas.bodies.size(), 32));
    m_update.set("uBodyCount", bodyCount);
    if (bodyCount > 0) {
        std::vector<glm::vec4> bodies(gas.bodies.begin(), gas.bodies.begin() + bodyCount);
        for (glm::vec4& b : bodies) b *= static_cast<float>(k);  // position and mass
        glProgramUniform4fv(m_update.id(), m_update.location("uBodies"), bodyCount, &bodies[0].x);
    }
    m_update.set("uLightDelay", delayed ? 1 : 0);
    m_update.set("uCameraPos", gas.cameraPosition * static_cast<float>(k));
    m_update.set("uCameraR", static_cast<float>(gas.cameraR * k));
    m_update.set("uTime", static_cast<float>(gas.time * k));
    m_update.set("uImpactDims", glm::ivec2(kImpactTheta, kImpactLongitude));
    m_update.set("uSurfaceShift", static_cast<float>(gas.surfaceShift));
    m_update.set("uSurfacePhase", static_cast<float>(gas.surfacePhase));
    // Last frame's gas, read while this frame's is accumulated.
    glBindTextureUnit(0, m_density);
    glBindTextureUnit(1, m_velocity);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_particles);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_grid);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, m_stats);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, m_impacts);
    glDispatchCompute((m_used + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    // Fold this frame's impacts into the fading hot-spot map.
    if (gas.surface > 0.0 && m_impactResolve.valid()) {
        m_impactResolve.use();
        m_impactResolve.set("uImpactDims", glm::ivec2(kImpactTheta, kImpactLongitude));
        m_impactResolve.set("uDecay", static_cast<float>(std::exp(-std::max(dt, 0.0) / std::max(gas.impactFadeTime, 1e-6))));
        m_impactResolve.set("uGain", 0.02f);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, m_impacts);
        glBindImageTexture(0, m_impactMap, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32F);
        glDispatchCompute((kImpactTheta + 7) / 8, (kImpactLongitude + 7) / 8, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_TEXTURE_FETCH_BARRIER_BIT);
    }

    // The gas as it is now, for the physics.
    resolve(m_density, m_velocity, GL_RGBA32F, true);

    // The gas as the camera sees it.
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glClearNamedBufferData(m_grid, GL_R32I, GL_RED_INTEGER, GL_INT, nullptr);
    m_deposit.use();
    m_deposit.set("uCount", m_used);
    m_deposit.set("uSpin", static_cast<float>(spin * k));
    m_deposit.set("uMass", static_cast<float>(k));
    m_deposit.set("uSeen", delayed ? 1 : 0);
    m_deposit.set("uGridDims", gridDims());
    m_deposit.set("uGridLogR", gridLogR());
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_particles);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_grid);
    glDispatchCompute((m_used + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    resolve(m_drawDensity, m_drawVelocity, GL_RGBA16F, false);
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_TEXTURE_UPDATE_BARRIER_BIT);

    // The alive count is for display only; reading it back stalls briefly, so do it rarely.
    if (++m_frame % 20 == 0) {
        int stats[2] = {0, 0};
        glGetNamedBufferSubData(m_stats, 0, sizeof(stats), stats);
        m_alive = stats[0];
        m_accreted = stats[1] / 256.0 / kObjectMass;
        if (m_alive == 0) clear();
    }
}

std::string ParticleSystem::describe() const {
    if (!active()) return "no debris";
    std::vector<GpuParticle> data(static_cast<size_t>(m_used));
    glGetNamedBufferSubData(m_particles, 0, sizeof(GpuParticle) * data.size(), data.data());
    std::vector<float> rs, ts, seen;
    int gone = 0;
    for (const GpuParticle& p : data) {
        if (p.state.w <= 0.0f) continue;
        if (p.state.w < 1.0f) ++gone;
        rs.push_back(static_cast<float>(p.state.x * m_scale));  // back to M
        ts.push_back(p.velocity.z);
        if (p.seen.w > p.seenVelocity.z && p.seen.w < p.seenVelocity.w) seen.push_back(static_cast<float>(p.seen.x * m_scale));
    }
    if (rs.empty()) return "no live debris";
    auto quantile = [](std::vector<float> v, double q) {
        std::sort(v.begin(), v.end());
        return v[static_cast<size_t>(q * (v.size() - 1))];
    };
    char text[384];
    int n = std::snprintf(text, sizeof(text),
                          "%zu parcels (%d gone, light in flight): r 10/50/90%% = %.1f/%.1f/%.1f M, "
                          "T 10/50/90/max = %.0f/%.0f/%.0f/%.0f K",
                          rs.size(), gone, quantile(rs, 0.1), quantile(rs, 0.5), quantile(rs, 0.9), quantile(ts, 0.1),
                          quantile(ts, 0.5), quantile(ts, 0.9), quantile(ts, 1.0));
    if (!seen.empty() && n > 0) {
        std::snprintf(text + n, sizeof(text) - n, "; seen: %zu, r 10/50/90%% = %.1f/%.1f/%.1f M", seen.size(),
                      quantile(seen, 0.1), quantile(seen, 0.5), quantile(seen, 0.9));
    }
    return text;
}

void ParticleSystem::hotReload() {
    m_update.reloadIfChanged();
    m_deposit.reloadIfChanged();
    m_resolve.reloadIfChanged();
    m_respin.reloadIfChanged();
    m_impactResolve.reloadIfChanged();
}

std::string ParticleSystem::errors() const {
    return m_update.error() + m_deposit.error() + m_resolve.error() + m_respin.error() + m_impactResolve.error();
}
