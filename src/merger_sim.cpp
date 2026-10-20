#include "merger_sim.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdint>

using namespace merger;

namespace {
constexpr double kEta = 1.2;            // smoothing length / particle spacing: ~58 neighbours
constexpr double kSofteningFactor = 1.2;  // gravitational softening, in core smoothing lengths
constexpr int kRelaxSteps = 300;
}  // namespace



bool MergerSim::init(const std::filesystem::path& shaderDir, std::string& error) {
    bool ok = m_grid.buildCompute(shaderDir, "merger_grid.comp");
    ok &= m_density.buildCompute(shaderDir, "merger_density.comp");
    ok &= m_force.buildCompute(shaderDir, "merger_force.comp");
    ok &= m_integrate.buildCompute(shaderDir, "merger_integrate.comp");
    ok &= m_reduce.buildCompute(shaderDir, "merger_reduce.comp");
    ok &= m_deposit.buildCompute(shaderDir, "merger_deposit.comp");
    ok &= m_resolve.buildCompute(shaderDir, "merger_resolve.comp");
    if (!ok) error = errors();

    glCreateBuffers(1, &m_head);
    glNamedBufferStorage(m_head, sizeof(int) * kHashSize, nullptr, GL_DYNAMIC_STORAGE_BIT);
    glCreateBuffers(1, &m_state);
    glNamedBufferStorage(m_state, sizeof(GpuState), nullptr, GL_DYNAMIC_STORAGE_BIT);
    glCreateBuffers(1, &m_history);
    glNamedBufferStorage(m_history, sizeof(GpuSample) * kHistorySize, nullptr, GL_DYNAMIC_STORAGE_BIT);
    const size_t cells = std::max<size_t>(static_cast<size_t>(kVolumeDims) * kVolumeDims * kVolumeDims * 8,
                                          static_cast<size_t>(kMapDims) * kMapDims * 2);
    glCreateBuffers(1, &m_cells);
    glNamedBufferStorage(m_cells, sizeof(int) * cells, nullptr, GL_DYNAMIC_STORAGE_BIT);
    glCreateBuffers(1, &m_counts);
    glNamedBufferStorage(m_counts, sizeof(int) * kCountDims * kCountDims * kCountDims, nullptr, GL_DYNAMIC_STORAGE_BIT);

    auto volume = [](GLuint& tex) {
        glCreateTextures(GL_TEXTURE_3D, 1, &tex);
        glTextureStorage3D(tex, 1, GL_RGBA16F, kVolumeDims, kVolumeDims, kVolumeDims);
        glTextureParameteri(tex, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTextureParameteri(tex, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        for (GLenum w : {GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T, GL_TEXTURE_WRAP_R})
            glTextureParameteri(tex, w, GL_CLAMP_TO_EDGE);
    };
    volume(m_volumeA);
    volume(m_volumeB);
    glCreateTextures(GL_TEXTURE_2D, 1, &m_map);
    glTextureStorage2D(m_map, 1, GL_RG32F, kMapDims, kMapDims);
    glTextureParameteri(m_map, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(m_map, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(m_map, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(m_map, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glCreateQueries(GL_TIME_ELAPSED, 1, &m_query);
    return ok;
}

void MergerSim::shutdown() {
    const GLuint buffers[] = {m_particles, m_head, m_next, m_state, m_aux, m_history, m_cells, m_counts};
    for (GLuint b : buffers)
        if (b) glDeleteBuffers(1, &b);
    const GLuint textures[] = {m_volumeA, m_volumeB, m_map};
    for (GLuint t : textures)
        if (t) glDeleteTextures(1, &t);
    if (m_query) glDeleteQueries(1, &m_query);
    m_particles = m_head = m_next = m_state = m_aux = m_history = m_cells = m_counts = 0;
    m_volumeA = m_volumeB = m_map = m_query = 0;
    m_count = 0;
}

void MergerSim::hotReload() {
    for (ShaderProgram* p : {&m_grid, &m_density, &m_force, &m_integrate, &m_reduce, &m_deposit, &m_resolve})
        p->reloadIfChanged();
}

std::string MergerSim::errors() const {
    std::string e;
    for (const ShaderProgram* p : {&m_grid, &m_density, &m_force, &m_integrate, &m_reduce, &m_deposit, &m_resolve})
        e += p->error();
    return e;
}

void MergerSim::bindAll() {
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, m_particles);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_head);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_next);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, m_state);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, m_aux);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, m_history);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, m_cells);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, m_counts);
}

void MergerSim::setCommon(ShaderProgram& p) {
    p.set("uCount", m_count);
    p.set("uCellSize", static_cast<float>(m_cell));
    p.set("uHashMask", kHashSize - 1);
}

// --- Building the stars ---------------------------------------------------------------------

void MergerSim::reset(const Params& params) {
    m_params = params;
    m_radius = params.radiusKm / kKm;
    m_K = polytropeK(m_radius);
    const double heavy = std::max(params.m1, params.m2);
    m_massRef = heavy / params.perStar;
    // Equal particle masses in both stars (SPH is happiest so).
    const int n0 = std::max(1000, static_cast<int>(std::lround(params.m1 / m_massRef)));
    const int n1 = std::max(1000, static_cast<int>(std::lround(params.m2 / m_massRef)));
    const double pm0 = params.m1 / n0, pm1 = params.m2 / n1;
    m_count = n0 + n1;
    m_hCore = kEta * std::cbrt(m_massRef / polytropeCentralDensity(heavy, m_radius));
    m_cell = 2.0 * m_hCore;

    if (m_count > m_capacity) {
        for (GLuint* b : {&m_particles, &m_next, &m_aux})
            if (*b) glDeleteBuffers(1, b);
        m_capacity = m_count;
        glCreateBuffers(1, &m_particles);
        glNamedBufferStorage(m_particles, sizeof(GpuParticle) * m_capacity, nullptr, GL_DYNAMIC_STORAGE_BIT);
        glCreateBuffers(1, &m_next);
        glNamedBufferStorage(m_next, sizeof(int) * m_capacity, nullptr, GL_DYNAMIC_STORAGE_BIT);
        glCreateBuffers(1, &m_aux);
        glNamedBufferStorage(m_aux, sizeof(glm::vec4) * m_capacity, nullptr, GL_DYNAMIC_STORAGE_BIT);
    }
    glClearNamedBufferData(m_aux, GL_RGBA32F, GL_RGBA, GL_FLOAT, nullptr);

    // Far apart (so they hardly pull on each other) while they settle.
    std::vector<GpuParticle> ps;
    ps.reserve(m_count);
    const double apart = 40.0 * m_radius;
    for (int star = 0; star < 2; ++star) {
        const int n = star == 0 ? n0 : n1;
        const double m = star == 0 ? params.m1 : params.m2;
        const double pm = star == 0 ? pm0 : pm1;
        const glm::dvec3 offset(star == 0 ? -apart : apart, 0.0, 0.0);
        for (const glm::dvec3& p : starParticles(n, m_radius)) {
            GpuParticle g{};
            g.pos = glm::vec4(glm::vec3(p + offset), static_cast<float>(pm));
            const double rho = std::max(polytropeDensity(m, m_radius, glm::length(p)), 0.02 * polytropeCentralDensity(m, m_radius));
            g.vel = glm::vec4(0.0f, 0.0f, 0.0f, static_cast<float>(kEta * std::cbrt(pm / rho)));
            g.thermo = glm::vec4(0.0f, static_cast<float>(rho), 0.0f, static_cast<float>(star));
            g.extra = glm::vec4(glm::vec3(p / m_radius), 0.0f);
            ps.push_back(g);
        }
    }
    glNamedBufferSubData(m_particles, 0, sizeof(GpuParticle) * ps.size(), ps.data());
    GpuState st{};
    glNamedBufferSubData(m_state, 0, sizeof(GpuState), &st);
    m_historyRead = 0;

    // Settle each star into equilibrium: damped steps, staying cold.
    bindAll();
    computeForces();
    for (int k = 0; k < kRelaxSteps; ++k) {
        m_reduce.use();
        setCommon(m_reduce);
        m_reduce.set("uDtMax", static_cast<float>(0.5 * m_hCore));
        m_reduce.set("uMergeSeparation", 0.0f);
        m_reduce.set("uHistorySize", kHistorySize);
        m_reduce.set("uRecord", 0);
        glDispatchCompute(1, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        m_integrate.use();
        setCommon(m_integrate);
        m_integrate.set("uStage", 0);
        glDispatchCompute((m_count + 255) / 256, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        computeForces();
        m_integrate.use();
        m_integrate.set("uStage", 3);
        m_integrate.set("uDamping", 0.9f);
        glDispatchCompute((m_count + 255) / 256, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }

    // Each parcel's place in its settled star becomes its label (it travels with the
    // matter: the cinematic view textures the stars with it).
    ps = readParticles();
    glm::dvec3 com[2] = {glm::dvec3(0.0), glm::dvec3(0.0)};
    double mass[2] = {0.0, 0.0};
    for (const GpuParticle& p : ps) {
        const int s = p.thermo.w < 0.5f ? 0 : 1;
        com[s] += static_cast<double>(p.pos.w) * glm::dvec3(p.pos);
        mass[s] += p.pos.w;
    }
    for (int s = 0; s < 2; ++s) com[s] /= std::max(mass[s], 1e-30);
    const float uth = static_cast<float>(thermalEnergyFromTemperature(params.initialTemperature));
    for (GpuParticle& p : ps) {
        const int s = p.thermo.w < 0.5f ? 0 : 1;
        p.extra = glm::vec4(glm::vec3((glm::dvec3(p.pos) - com[s]) / m_radius), 0.0f);
        p.vel = glm::vec4(0.0f, 0.0f, 0.0f, p.vel.w);
        p.acc = glm::vec4(0.0f);
        p.thermo.x = uth;
    }
    glNamedBufferSubData(m_particles, 0, sizeof(GpuParticle) * ps.size(), ps.data());
    st = GpuState{};
    glNamedBufferSubData(m_state, 0, sizeof(GpuState), &st);

    // Peters' inspiral from the starting separation to where the hydrodynamics takes over.
    const double M = params.m1 + params.m2;
    const double a0 = params.startKm / kKm, aHydro = std::min(params.hydroKm / kKm, a0);
    const double T = coalescenceTime(a0, params.m1, params.m2);
    m_hydroStart = T - coalescenceTime(aHydro, params.m1, params.m2);
    m_orbitDt = std::max(m_hydroStart / 200000.0, 0.05);
    m_orbit.clear();
    double phase = 0.0;
    for (double t = 0.0;; t += m_orbitDt) {
        const double a = separationAfter(a0, t, params.m1, params.m2);
        m_orbit.push_back({a, phase});
        if (t >= m_hydroStart) break;
        // Midpoint rule for d phase / dt = sqrt(M / a^3).
        const double am = separationAfter(a0, t + 0.5 * m_orbitDt, params.m1, params.m2);
        phase += m_orbitDt * std::sqrt(M / (am * am * am));
    }
    m_mergeSeparation = 1.2 * m_radius;

    m_phase = Phase::Inspiral;
    m_time = 0.0;
    m_mergeTime = -1.0;
    m_coastStart = -1.0;
    m_holeTime = -1.0;
    m_hole = BlackHole{};
    m_steps = 0;
    m_sortedAt = -1000;
    m_lastDt = 0.3 * m_hCore;
    m_waves.clear();
    m_gwFrequency = 0.0;
    m_unbound = 0.0;
    m_aliveMass = M;
    m_ejectaSky.assign(32 * 16, 0.0f);
    m_maxDensity = polytropeCentralDensity(heavy, m_radius) * kGramsPerCc;
    m_maxTemperature = params.initialTemperature;
    carry(0.0);
}

std::vector<MergerSim::GpuParticle> MergerSim::readParticles() const {
    std::vector<GpuParticle> ps(static_cast<size_t>(m_count));
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
    glGetNamedBufferSubData(m_particles, 0, sizeof(GpuParticle) * ps.size(), ps.data());
    return ps;
}

// --- Stepping ------------------------------------------------------------------------------------

void MergerSim::computeForces() {
    m_grid.use();
    setCommon(m_grid);
    m_grid.set("uStage", 0);
    glDispatchCompute(kHashSize / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    m_grid.set("uStage", 1);
    glDispatchCompute((m_count + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    m_density.use();
    setCommon(m_density);
    m_density.set("uK", static_cast<float>(m_K));
    m_density.set("uGammaTh", static_cast<float>(kGammaThermal));
    m_density.set("uHMin", static_cast<float>(0.2 * m_hCore));
    m_density.set("uHMax", static_cast<float>(1.5 * m_cell));
    m_density.set("uEta", static_cast<float>(kEta));
    glDispatchCompute((m_count + 127) / 128, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    m_force.use();
    setCommon(m_force);
    m_force.set("uAlpha", 1.0f);
    m_force.set("uBeta", 2.0f);
    m_force.set("uSoftening", static_cast<float>(kSofteningFactor * m_hCore));
    m_force.set("uGammaTh", static_cast<float>(kGammaThermal));
    m_force.set("uGravity", 1);
    glDispatchCompute((m_count + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

void MergerSim::dispatchStep(double dtMax, bool record) {
    m_reduce.use();
    setCommon(m_reduce);
    m_reduce.set("uDtMax", static_cast<float>(dtMax));
    m_reduce.set("uMergeSeparation", static_cast<float>(m_mergeSeparation));
    m_reduce.set("uHistorySize", kHistorySize);
    m_reduce.set("uRecord", record ? 1 : 0);
    glDispatchCompute(1, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    m_integrate.use();
    setCommon(m_integrate);
    m_integrate.set("uStage", 0);
    glDispatchCompute((m_count + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    computeForces();

    m_integrate.use();
    m_integrate.set("uStage", 1);
    glDispatchCompute((m_count + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

void MergerSim::carry(double t) {
    const double M = m_params.m1 + m_params.m2;
    const double u = std::clamp(t / m_orbitDt, 0.0, static_cast<double>(m_orbit.size() - 1));
    const size_t k = std::min(static_cast<size_t>(u), m_orbit.size() - 1);
    const size_t k1 = std::min(k + 1, m_orbit.size() - 1);
    const double f = u - static_cast<double>(k);
    const double a = m_orbit[k].x * (1.0 - f) + m_orbit[k1].x * f;
    const double phi = m_orbit[k].y * (1.0 - f) + m_orbit[k1].y * f;
    const double omega = std::sqrt(M / (a * a * a));
    const double adot = -64.0 / 5.0 * m_params.m1 * m_params.m2 * M / (a * a * a);
    const glm::dvec3 n(std::cos(phi), std::sin(phi), 0.0), tangent(-n.y, n.x, 0.0);
    const glm::dvec3 rel = a * n, relV = adot * n + a * omega * tangent;
    m_center[0] = m_params.m2 / M * rel;
    m_center[1] = -m_params.m1 / M * rel;
    m_velocity[0] = m_params.m2 / M * relV;
    m_velocity[1] = -m_params.m1 / M * relV;
    m_separation = a;

    bindAll();
    m_integrate.use();
    setCommon(m_integrate);
    m_integrate.set("uStage", 2);
    m_integrate.set("uRadius", static_cast<float>(m_radius));
    for (int s = 0; s < 2; ++s) {
        const std::string c = "uCenter[" + std::to_string(s) + "]", v = "uVelocity[" + std::to_string(s) + "]";
        m_integrate.set(c.c_str(), glm::vec3(m_center[s]));
        m_integrate.set(v.c_str(), glm::vec3(m_velocity[s]));
    }
    glDispatchCompute((m_count + 255) / 256, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

void MergerSim::startHydro() {
    carry(m_hydroStart);
    m_time = m_hydroStart;
    GpuState st{};
    st.time = m_time;
    glNamedBufferSubData(m_state, 0, sizeof(GpuState), &st);
    bindAll();
    computeForces();
    m_phase = Phase::Hydro;
}

void MergerSim::readState() {
    GpuState st{};
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
    glGetNamedBufferSubData(m_state, 0, sizeof(GpuState), &st);
    m_time = st.time;
    m_lastDt = std::max(static_cast<double>(st.step.y), 1e-3);
    m_steps = static_cast<int>(st.step.z);
    m_center[0] = glm::dvec3(st.star[0]);
    m_center[1] = glm::dvec3(st.star[2]);
    m_velocity[0] = glm::dvec3(st.star[1]);
    m_velocity[1] = glm::dvec3(st.star[3]);
    if (st.star[0].w > 0.0f && st.star[2].w > 0.0f) m_separation = glm::length(m_center[0] - m_center[1]);
    m_aliveMass = static_cast<double>(st.star[0].w) + st.star[2].w;
    m_maxDensity = st.iddot1.z * kGramsPerCc;
    m_maxTemperature = std::pow(10.0, static_cast<double>(st.iddot1.w));
    m_rms = st.extent.x;
    m_extent = st.extent.y;
    m_unbound = st.extent.z;
    if (st.step.w > 0.5f && m_mergeTime < 0.0) m_mergeTime = m_time;

    // New wave samples.
    const int written = st.counters.x;
    const int first = std::max(m_historyRead, written - kHistorySize);
    if (written > first) {
        std::vector<GpuSample> samples(static_cast<size_t>(written - first));
        // The ring may wrap: read in up to two pieces.
        int k = first;
        size_t out = 0;
        while (k < written) {
            const int slot = k % kHistorySize;
            const int n = std::min(written - k, kHistorySize - slot);
            glGetNamedBufferSubData(m_history, sizeof(GpuSample) * slot, sizeof(GpuSample) * n, &samples[out]);
            out += static_cast<size_t>(n);
            k += n;
        }
        for (const GpuSample& s : samples) {
            m_waves.push_back({s.t, {s.a.x, s.a.y, s.a.z, s.a.w, s.b.x, s.b.y}});
        }
        m_historyRead = written;
    }
    // The wave's frequency from how fast its phase (seen from above) turns.
    if (m_waves.size() > 8) {
        const WavePoint& p = m_waves.back();
        const WavePoint& q = m_waves[m_waves.size() - 8];
        const double ph1 = std::atan2(p.iddot[3], 0.5 * (p.iddot[0] - p.iddot[1]));
        const double ph0 = std::atan2(q.iddot[3], 0.5 * (q.iddot[0] - q.iddot[1]));
        const double d = std::remainder(ph1 - ph0, 2.0 * kPi);
        if (p.t > q.t) m_gwFrequency = std::abs(d) / (2.0 * kPi * (p.t - q.t));
    }
}

// Neighbours far apart in memory make the neighbour search slow (every candidate is a
// cache miss): every so often the particles are put in Morton order (a Z-curve through the
// grid cells), so particles near each other in space sit near each other in memory too.
void MergerSim::sortParticles() {
    std::vector<GpuParticle> ps = readParticles();
    auto spread = [](uint64_t v) {
        v &= 0x1fffff;
        v = (v | v << 32) & 0x1f00000000ffffULL;
        v = (v | v << 16) & 0x1f0000ff0000ffULL;
        v = (v | v << 8) & 0x100f00f00f00f00fULL;
        v = (v | v << 4) & 0x10c30c30c30c30c3ULL;
        v = (v | v << 2) & 0x1249249249249249ULL;
        return v;
    };
    std::vector<std::pair<uint64_t, int>> keys(ps.size());
    for (size_t i = 0; i < ps.size(); ++i) {
        uint64_t code = ~0ULL;  // swallowed parcels go last
        if (ps[i].thermo.w >= 0.0f) {
            const glm::dvec3 c = glm::floor(glm::dvec3(ps[i].pos) / m_cell) + 1048576.0;
            const glm::dvec3 k = glm::clamp(c, glm::dvec3(0.0), glm::dvec3(2097151.0));
            code = spread(static_cast<uint64_t>(k.x)) | spread(static_cast<uint64_t>(k.y)) << 1 |
                   spread(static_cast<uint64_t>(k.z)) << 2;
        }
        keys[i] = {code, static_cast<int>(i)};
    }
    std::sort(keys.begin(), keys.end());
    std::vector<GpuParticle> sorted(ps.size());
    for (size_t i = 0; i < ps.size(); ++i) sorted[i] = ps[static_cast<size_t>(keys[i].second)];
    glNamedBufferSubData(m_particles, 0, sizeof(GpuParticle) * sorted.size(), sorted.data());
    m_sortedAt = m_steps;
}

void MergerSim::advance(double target, int maxSteps) {
    if (!ready()) return;
    bindAll();
    if (m_phase == Phase::Inspiral) {
        const double M = m_params.m1 + m_params.m2;
        if (target < m_hydroStart) {
            m_time = target;
            carry(target);
            // The waves of two point masses on a circular orbit, as they were then.
            const double u = std::clamp(target / m_orbitDt, 0.0, static_cast<double>(m_orbit.size() - 1));
            const size_t k = static_cast<size_t>(u);
            const double a = m_orbit[k].x;
            m_gwFrequency = std::sqrt(M / (a * a * a)) / kPi;
            return;
        }
        startHydro();
    }
    if (m_phase == Phase::Coasting) {
        m_time = std::max(m_time, target);
        return;
    }
    glBeginQuery(GL_TIME_ELAPSED, m_query);
    int done = 0;
    while (done < maxSteps && m_time < target - 1e-9) {
        const double left = target - m_time;
        if (m_steps - m_sortedAt >= 16) {
            sortParticles();
            bindAll();
            computeForces();  // the per-particle rates (aux) follow the new order
        }
        const int n = std::clamp(static_cast<int>(std::ceil(left / (0.95 * m_lastDt))), 1,
                                 std::min(maxSteps - done, 16));
        const double dtMax = left / n;
        for (int k = 0; k < n; ++k) dispatchStep(dtMax, true);
        done += n;
        readState();
        if (m_phase == Phase::Hydro && merged() && m_time >= collapseTime()) {
            doCollapse();
            break;
        }
        if (m_phase == Phase::Collapsed && m_time >= m_holeTime + m_params.coastMs * 1e-3 / kSecond) {
            measureEjecta();
            m_phase = Phase::Coasting;
            m_coastStart = m_time;
            break;
        }
    }
    glEndQuery(GL_TIME_ELAPSED);
    GLuint64 ns = 0;
    glGetQueryObjectui64v(m_query, GL_QUERY_RESULT, &ns);
    m_gpuMs = static_cast<float>(ns * 1e-6);
}

double MergerSim::sinceMerger() const {
    if (m_mergeTime >= 0.0) return m_time - m_mergeTime;
    // Predicted from Peters' law: when the centres would be m_mergeSeparation apart.
    const double a = std::max(m_separation, m_mergeSeparation);
    return -(coalescenceTime(a, m_params.m1, m_params.m2) -
             coalescenceTime(m_mergeSeparation, m_params.m1, m_params.m2));
}

double MergerSim::holeMass() const {
    return hasHole() ? totalMass() - m_aliveMass : 0.0;
}

// --- Collapse and ejecta -----------------------------------------------------------------------

void MergerSim::doCollapse() {
    std::vector<GpuParticle> ps = readParticles();
    const size_t n = ps.size();
    std::vector<glm::dvec3> x(n), v(n);
    std::vector<double> m(n), phi(n);
    std::vector<bool> alive(n), swallowed;
    double rhoMax = 0.0;
    for (const GpuParticle& p : ps)
        if (p.thermo.w >= 0.0f) rhoMax = std::max(rhoMax, static_cast<double>(p.thermo.y));
    glm::dvec3 center(0.0);
    double cm = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const GpuParticle& p = ps[i];
        x[i] = glm::dvec3(p.pos);
        v[i] = glm::dvec3(p.vel);
        m[i] = p.pos.w;
        phi[i] = p.extra.w;
        alive[i] = p.thermo.w >= 0.0f;
        if (alive[i] && p.thermo.y > 0.3 * rhoMax) {
            center += m[i] * x[i];
            cm += m[i];
        }
    }
    center /= std::max(cm, 1e-30);
    m_hole = collapse(x, v, m, phi, alive, center, swallowed);
    for (size_t i = 0; i < n; ++i)
        if (swallowed[i]) ps[i].thermo.w = -1.0f;
    glNamedBufferSubData(m_particles, 0, sizeof(GpuParticle) * n, ps.data());
    const glm::vec4 hole(glm::vec3(m_hole.position), static_cast<float>(m_hole.mass));
    const glm::vec4 holeV(glm::vec3(m_hole.velocity), 0.0f);
    glNamedBufferSubData(m_state, offsetof(GpuState, hole), sizeof(glm::vec4), &hole);
    glNamedBufferSubData(m_state, offsetof(GpuState, holeV), sizeof(glm::vec4), &holeV);
    m_holeTime = m_time;
    m_ring = ringdown(m_hole.mass, m_hole.spin);
    if (!m_waves.empty()) {
        const auto& I = m_waves.back().iddot;
        m_ringAmplitude[0] = std::hypot(0.5 * (I[0] - I[1]), I[3]);
        m_ringAmplitude[1] = std::atan2(I[3], 0.5 * (I[0] - I[1]));
    }
    m_phase = Phase::Collapsed;
    std::printf("Merger: the remnant collapsed %.2f ms after the merger into a black hole of %.3f Suns, spin %.2f; "
                "torus %.3f Suns, unbound %.4f Suns\n",
                (m_time - m_mergeTime) * kSecond * 1e3, m_hole.mass, m_hole.spin, m_hole.diskMass, m_hole.unboundMass);
    bindAll();
    computeForces();
}

void MergerSim::measureEjecta() {
    std::vector<GpuParticle> ps = readParticles();
    m_ejectaSky.assign(32 * 16, 0.0f);
    double mass = 0.0, speed = 0.0;
    for (const GpuParticle& p : ps) {
        if (p.thermo.w < 0.0f) continue;
        const glm::dvec3 vel(p.vel), pos(p.pos);
        const double e = 0.5 * glm::dot(vel, vel) + p.extra.w;
        if (e <= 0.0 || glm::dot(vel, pos) <= 0.0) continue;
        const double s = glm::length(vel);
        const double th = std::acos(std::clamp(vel.z / s, -1.0, 1.0));
        const double ph = std::atan2(vel.y, vel.x) + kPi;
        const int i = std::clamp(static_cast<int>(ph / (2.0 * kPi) * 32.0), 0, 31);
        const int j = std::clamp(static_cast<int>(th / kPi * 16.0), 0, 15);
        m_ejectaSky[j * 32 + i] += p.pos.w;
        mass += p.pos.w;
        speed += p.pos.w * std::sqrt(2.0 * e);  // speed at infinity
    }
    if (mass > 0.0) {
        for (float& f : m_ejectaSky) f = static_cast<float>(f / mass);
        m_ejectaSpeed = speed / mass;
    }
    std::printf("Merger: %.4f Suns unbound, coasting out at %.3f c on average\n", mass, m_ejectaSpeed);
}

// --- Waves -------------------------------------------------------------------------------------------

Strain MergerSim::strainAt(double t, const glm::dvec3& n, double distance) const {
    return strain(iddotAt(t), n, distance);
}

std::array<double, 6> MergerSim::iddotAt(double t) const {
    std::array<double, 6> I{};
    if (m_orbit.empty()) return I;
    const double M = m_params.m1 + m_params.m2, mu = m_params.m1 * m_params.m2 / M;
    if (t < m_hydroStart || m_waves.empty() || t < m_waves.front().t) {
        // Two point masses on Peters' inspiral.
        const double u = std::clamp(t / m_orbitDt, 0.0, static_cast<double>(m_orbit.size() - 1));
        const size_t k = std::min(static_cast<size_t>(u), m_orbit.size() - 1);
        const size_t k1 = std::min(k + 1, m_orbit.size() - 1);
        const double f = u - static_cast<double>(k);
        const double a = m_orbit[k].x * (1.0 - f) + m_orbit[k1].x * f;
        const double phi = m_orbit[k].y * (1.0 - f) + m_orbit[k1].y * f;
        const double w2 = M / (a * a * a);
        const double A = 2.0 * mu * a * a * w2;
        I = {-A * std::cos(2 * phi), A * std::cos(2 * phi), 0.0, -A * std::sin(2 * phi), 0.0, 0.0};
    } else {
        auto it = std::lower_bound(m_waves.begin(), m_waves.end(), t,
                                   [](const WavePoint& p, double x) { return p.t < x; });
        if (it == m_waves.end()) {
            I = m_waves.back().iddot;
        } else if (it == m_waves.begin()) {
            I = it->iddot;
        } else {
            const WavePoint& b = *it;
            const WavePoint& a = *(it - 1);
            const double f = (t - a.t) / std::max(b.t - a.t, 1e-12);
            for (int c = 0; c < 6; ++c) I[c] = a.iddot[c] * (1.0 - f) + b.iddot[c] * f;
        }
    }
    if (m_holeTime >= 0.0 && t > m_holeTime) {
        // The new hole rings down from the remnant's last quadrupole.
        const double dt = t - m_holeTime;
        const double A = m_ringAmplitude[0] * std::exp(-dt / m_ring.tau);
        const double ph = m_ringAmplitude[1] + m_ring.omega * dt;
        I[0] += A * std::cos(ph);
        I[1] -= A * std::cos(ph);
        I[3] += A * std::sin(ph);
    }
    return I;
}

// --- Drawing grids -------------------------------------------------------------------------------

void MergerSim::deposit(const glm::vec3& center, float half, const glm::vec2& mapCenter, float mapHalf) {
    if (!ready()) return;
    bindAll();
    const float coast = m_phase == Phase::Coasting ? static_cast<float>(m_time - m_coastStart) : 0.0f;
    m_boxMin = center - glm::vec3(half);
    m_boxSize = 2.0f * half;
    m_mapMin = mapCenter - glm::vec2(mapHalf);
    m_mapSize = 2.0f * mapHalf;
    const float boxCell = m_boxSize / kVolumeDims, mapCell = m_mapSize / kMapDims;
    // Mass / length^3 and mass / length^2 in g/cm^3 and g/cm^2.
    const float densityUnit = static_cast<float>(kGramsPerCc);
    const float columnUnit = static_cast<float>(kSunGrams / std::pow(kKm * 1e5, 2));

    // First count the particles on a coarse grid over everything drawn (how sparse the gas is).
    const float countHalf = std::max(half, mapHalf);
    const float countCell = 2.0f * countHalf / kCountDims;
    auto setDeposit = [&](int stage) {
        m_deposit.use();
        setCommon(m_deposit);
        m_deposit.set("uStage", stage);
        m_deposit.set("uCountMin", glm::vec3(-countHalf));
        m_deposit.set("uCountCell", countCell);
        m_deposit.set("uCountDims", kCountDims);
        m_deposit.set("uMassRef", static_cast<float>(m_massRef));
        m_deposit.set("uCoast", coast);
    };
    glClearNamedBufferData(m_counts, GL_R32I, GL_RED_INTEGER, GL_INT, nullptr);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    setDeposit(2);
    glDispatchCompute((m_count + 127) / 128, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    for (int stage = 0; stage < 2; ++stage) {
        glClearNamedBufferData(m_cells, GL_R32I, GL_RED_INTEGER, GL_INT, nullptr);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        setDeposit(stage);
        m_deposit.set("uBoxMin", m_boxMin);
        m_deposit.set("uBoxCell", boxCell);
        m_deposit.set("uBoxDims", kVolumeDims);
        m_deposit.set("uMapMin", m_mapMin);
        m_deposit.set("uMapCell", mapCell);
        m_deposit.set("uMapDims", kMapDims);
        glDispatchCompute((m_count + 127) / 128, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        m_resolve.use();
        m_resolve.set("uStage", stage);
        m_resolve.set("uBoxDims", kVolumeDims);
        m_resolve.set("uBoxCell", boxCell);
        m_resolve.set("uMapDims", kMapDims);
        m_resolve.set("uMapCell", mapCell);
        m_resolve.set("uMassRef", static_cast<float>(m_massRef));
        m_resolve.set("uDensityUnit", densityUnit);
        m_resolve.set("uColumnUnit", columnUnit);
        glBindImageTexture(0, m_volumeA, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA16F);
        glBindImageTexture(1, m_volumeB, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA16F);
        glBindImageTexture(2, m_map, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
        if (stage == 0) glDispatchCompute(kVolumeDims / 4, kVolumeDims / 4, kVolumeDims / 4);
        else glDispatchCompute(kMapDims / 4, kMapDims / 4, 1);
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
}
