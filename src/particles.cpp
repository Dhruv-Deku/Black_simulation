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
    glNamedBufferStorage(m_stats, sizeof(int), nullptr, GL_DYNAMIC_STORAGE_BIT);

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
    GLuint buffers[] = {m_particles, m_grid, m_stats};
    glDeleteBuffers(3, buffers);
    GLuint textures[] = {m_density, m_velocity, m_drawDensity, m_drawVelocity, m_occupancy};
    glDeleteTextures(5, textures);
    m_particles = m_grid = m_stats = 0;
    m_density = m_velocity = m_drawDensity = m_drawVelocity = m_occupancy = 0;
    m_used = m_next = m_alive = 0;
}

void ParticleSystem::spawn(const std::vector<ParticleInit>& particles, double time, bool existedBefore) {
    if (particles.empty()) return;
    allocate();
    const float weight = kObjectMass / static_cast<float>(particles.size());
    const float now = static_cast<float>(time);
    const float created = existedBefore ? -1e30f : now;
    std::vector<GpuParticle> data;
    data.reserve(particles.size());
    for (const ParticleInit& p : particles) {
        data.push_back({glm::vec4(p.r, p.theta, p.phi, 1.0f), glm::vec4(p.vr, p.vtheta, p.temperature, weight),
                        glm::vec4(p.E, p.L, p.Q, 0.0f), glm::vec4(p.r, p.theta, p.phi, now),
                        glm::vec4(p.vr, p.vtheta, created, 1e30f)});
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

void ParticleSystem::respin(double oldSpin, double newSpin) {
    if (!active() || !m_respin.valid()) return;
    m_respin.use();
    m_respin.set("uCount", m_used);
    m_respin.set("uOldSpin", static_cast<float>(oldSpin));
    m_respin.set("uSpin", static_cast<float>(newSpin));
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
    glClearNamedBufferData(m_stats, GL_R32I, GL_RED_INTEGER, GL_INT, nullptr);
    const GLuint zero = 0;
    glClearTexImage(m_occupancy, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, &zero);

    m_update.use();
    m_update.set("uCount", m_used);
    m_update.set("uSpin", static_cast<float>(spin));
    m_update.set("uDt", static_cast<float>(dt));
    m_update.set("uGridDims", gridDims());
    m_update.set("uGridLogR", gridLogR());
    m_update.set("uGasPhysics", gas.physics ? 1 : 0);
    m_update.set("uCollisionRate", gas.collisionRate);
    m_update.set("uCoolingTime", std::max(gas.coolingTime, 1.0f));
    m_update.set("uFloorTemperature", gas.floorTemperature);
    const int bodyCount = static_cast<int>(std::min<size_t>(gas.bodies.size(), 32));
    m_update.set("uBodyCount", bodyCount);
    if (bodyCount > 0) {
        glProgramUniform4fv(m_update.id(), m_update.location("uBodies"), bodyCount, &gas.bodies[0].x);
    }
    m_update.set("uLightDelay", delayed ? 1 : 0);
    m_update.set("uCameraPos", gas.cameraPosition);
    m_update.set("uCameraR", gas.cameraR);
    m_update.set("uTime", static_cast<float>(gas.time));
    // Last frame's gas, read while this frame's is accumulated.
    glBindTextureUnit(0, m_density);
    glBindTextureUnit(1, m_velocity);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_particles);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_grid);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, m_stats);
    glDispatchCompute((m_used + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    // The gas as it is now, for the physics.
    resolve(m_density, m_velocity, GL_RGBA32F, true);

    // The gas as the camera sees it.
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glClearNamedBufferData(m_grid, GL_R32I, GL_RED_INTEGER, GL_INT, nullptr);
    m_deposit.use();
    m_deposit.set("uCount", m_used);
    m_deposit.set("uSpin", static_cast<float>(spin));
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
        glGetNamedBufferSubData(m_stats, 0, sizeof(int), &m_alive);
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
        rs.push_back(p.state.x);
        ts.push_back(p.velocity.z);
        if (p.seen.w > p.seenVelocity.z && p.seen.w < p.seenVelocity.w) seen.push_back(p.seen.x);
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
}

std::string ParticleSystem::errors() const {
    return m_update.error() + m_deposit.error() + m_resolve.error() + m_respin.error();
}
