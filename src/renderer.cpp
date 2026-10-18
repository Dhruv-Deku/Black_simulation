#include "renderer.h"
#include <cstdio>

#include "blackbody.h"
#include "kerr.h"

#include <algorithm>
#include <cmath>

#ifdef _MSC_VER
#pragma warning(push, 0)
#endif
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

namespace {
constexpr int kBlackbodySamples = 1024;
constexpr double kFlowPeriod = 300.0;  // M; how often the disk's turbulence pattern is renewed

// Cells binning path segments for the ray tracer, in log r from the debris grid's inner
// radius out to its outer one, or farther when orbits reach farther (moons and planets).
constexpr int kPathR = 64, kPathTheta = 48, kPathPhi = 96;
constexpr int kPathCells = kPathR * kPathTheta * kPathPhi;

struct GpuPathSegment {
    glm::vec4 a, b, color, dash;
};
static_assert(sizeof(GpuPathSegment) == 64, "std430 layout");

// The path cell containing a Cartesian point (render units), or -1 outside the grid
// (which reaches out to rMax).
int pathCell(const glm::vec3& p, double spin, double rMax) {
    auto [r, theta, phi] = kerr::fromCartesian(spin, p.x, p.y, p.z);
    const double lr = (std::log(r) - std::log(ParticleSystem::kGridRMin)) /
                      (std::log(rMax) - std::log(ParticleSystem::kGridRMin));
    if (lr < 0.0 || lr >= 1.0) return -1;
    double f = phi / (2.0 * kerr::kPi);
    f -= std::floor(f);
    const int x = std::min(static_cast<int>(lr * kPathR), kPathR - 1);
    const int y = std::clamp(static_cast<int>(theta / kerr::kPi * kPathTheta), 0, kPathTheta - 1);
    const int z = std::min(static_cast<int>(f * kPathPhi), kPathPhi - 1);
    return (z * kPathTheta + y) * kPathR + x;
}
}  // namespace

ShaderProgram& Renderer::traceProgram(int look, bool lenses, bool planets, bool outflow) {
    // look 0: a black hole; 1-6 (surface::Look): a star or planet at the center.
    look = std::clamp(look, 0, kTraceLooks - 1);
    const int index = ((look * 2 + (lenses ? 1 : 0)) * 2 + (planets ? 1 : 0)) * 2 + (outflow ? 1 : 0);
    ShaderProgram& program = m_trace[index];
    if (!program.valid() && program.error().empty()) {
        std::string defines;
        if (look > 0) defines += "#define CENTRAL_BODY 1\n#define SURFACE_LOOK " + std::to_string(look) + "\n";
        if (lenses) defines += "#define LENSES 1\n";
        if (planets) defines += "#define BODY_PLANETS 1\n";
        if (outflow) defines += "#define OUTFLOW 1\n";
        program.build(m_shaderDir, "fullscreen.vert", "trace.frag", defines);
    }
    return program;
}

bool Renderer::init(const std::filesystem::path& shaderDir, std::string& error) {
    m_shaderDir = shaderDir;
    bool ok = traceProgram(0, false, false, false).valid();
    ok &= m_bloomDown.build(shaderDir, "fullscreen.vert", "bloom_down.frag");
    ok &= m_bloomUp.build(shaderDir, "fullscreen.vert", "bloom_up.frag");
    ok &= m_tonemap.build(shaderDir, "fullscreen.vert", "tonemap.frag");
    // Keep going on failure: a fixed shader file is picked up by hotReload().
    if (!ok) error = shaderErrors();

    glCreateVertexArrays(1, &m_vao);

    glCreateBuffers(1, &m_bodyUbo);
    glNamedBufferStorage(m_bodyUbo, sizeof(GpuBodies), nullptr, GL_DYNAMIC_STORAGE_BIT);

    // 1D table of blackbody colors, indexed by log10(temperature).
    std::vector<float> table = blackbody::makeTable(kBlackbodySamples);
    glCreateTextures(GL_TEXTURE_2D, 1, &m_blackbodyTexture);
    glTextureStorage2D(m_blackbodyTexture, 1, GL_RGB32F, kBlackbodySamples, 1);
    glTextureSubImage2D(m_blackbodyTexture, 0, 0, 0, kBlackbodySamples, 1, GL_RGB, GL_FLOAT, table.data());
    glTextureParameteri(m_blackbodyTexture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(m_blackbodyTexture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(m_blackbodyTexture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(m_blackbodyTexture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Bound in place of the sky image until one is loaded.
    const unsigned char black[4] = {0, 0, 0, 255};
    glCreateTextures(GL_TEXTURE_2D, 1, &m_placeholderTexture);
    glTextureStorage2D(m_placeholderTexture, 1, GL_RGBA8, 1, 1);
    glTextureSubImage2D(m_placeholderTexture, 0, 0, 0, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, black);

    // Bound in place of the debris volumes while there is no debris.
    const float empty[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    glCreateTextures(GL_TEXTURE_3D, 1, &m_placeholderVolume);
    glTextureStorage3D(m_placeholderVolume, 1, GL_RGBA16F, 1, 1, 1);
    glTextureSubImage3D(m_placeholderVolume, 0, 0, 0, 0, 1, 1, 1, GL_RGBA, GL_FLOAT, empty);
    const GLuint nothing = 0;
    glCreateTextures(GL_TEXTURE_3D, 1, &m_placeholderOccupancy);
    glTextureStorage3D(m_placeholderOccupancy, 1, GL_R32UI, 1, 1, 1);
    glTextureSubImage3D(m_placeholderOccupancy, 0, 0, 0, 0, 1, 1, 1, GL_RED_INTEGER, GL_UNSIGNED_INT, &nothing);
    glTextureParameteri(m_placeholderOccupancy, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(m_placeholderOccupancy, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    // Past positions and 4-velocities of the bodies, for light-travel delay.
    glCreateTextures(GL_TEXTURE_2D, 1, &m_historyTexture);
    glTextureStorage2D(m_historyTexture, 1, GL_RGBA32F, History::kColumns, History::kRows * History::kSlots);
    glTextureParameteri(m_historyTexture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(m_historyTexture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glCreateQueries(GL_TIME_ELAPSED, kQueries, m_queries);

    // Path line buffers start out empty (every cell holds zero segments).
    glCreateBuffers(1, &m_pathSegments);
    glCreateBuffers(1, &m_pathCells);
    glCreateBuffers(1, &m_pathIndex);
    glCreateTextures(GL_TEXTURE_3D, 1, &m_pathOccupancy);
    glTextureStorage3D(m_pathOccupancy, 1, GL_R8UI, kPathR, kPathTheta, kPathPhi);
    glTextureParameteri(m_pathOccupancy, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(m_pathOccupancy, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    setPaths({}, 0.0, 1.0);
    loadPlanetMaps();
    m_sdoBake.build(shaderDir, "fullscreen.vert", "sdo_bake.frag");
    return ok;
}

void Renderer::setPaths(const std::vector<PathLine>& lines, double spin, double scale, bool farOrbits) {
    std::vector<GpuPathSegment> segments;
    std::vector<std::vector<uint32_t>> cells(kPathCells);
    const float inv = static_cast<float>(1.0 / scale);
    const double spinRender = spin / scale;
    double rMax = ParticleSystem::kGridRMax;
    for (const PathLine& line : lines) {
        if (!farOrbits) break;
        for (const glm::vec3& point : line.points) rMax = std::max(rMax, 1.05 * glm::length(point * inv));
    }
    m_pathRMax = rMax;
    for (const PathLine& line : lines) {
        const size_t n = line.points.size();
        if (n < 2) continue;
        float along = 0.0f;
        int prevCell = pathCell(line.points[0] * inv, spinRender, rMax);
        for (size_t i = 1; i < n; ++i) {
            const glm::vec3 a = line.points[i - 1] * inv;
            const glm::vec3 b = line.points[i] * inv;
            const float len = glm::length(b - a);
            const float fa = static_cast<float>(i - 1) / static_cast<float>(n - 1);
            const float fb = static_cast<float>(i) / static_cast<float>(n - 1);
            const int cell = pathCell(b, spinRender, rMax);
            const int mid = pathCell(0.5f * (a + b), spinRender, rMax);
            const uint32_t index = static_cast<uint32_t>(segments.size());
            segments.push_back({glm::vec4(a, glm::mix(line.startOpacity, line.endOpacity, fa)),
                                glm::vec4(b, glm::mix(line.startOpacity, line.endOpacity, fb)),
                                glm::vec4(line.color, line.dashed ? 1.0f : 0.0f),
                                glm::vec4(along, along + len, 0.0f, 0.0f)});
            along += len;
            for (int c : {prevCell, cell, mid}) {
                if (c >= 0 && (cells[c].empty() || cells[c].back() != index)) cells[c].push_back(index);
            }
            prevCell = cell;
        }
    }

    // Per cell: where its segments start in `indices`, and how many there are.
    std::vector<glm::uvec2> ranges(kPathCells);
    std::vector<uint32_t> indices;
    std::vector<uint8_t> occupied(kPathCells, 0);
    for (int c = 0; c < kPathCells; ++c) {
        ranges[c] = glm::uvec2(static_cast<uint32_t>(indices.size()), static_cast<uint32_t>(cells[c].size()));
        indices.insert(indices.end(), cells[c].begin(), cells[c].end());
        occupied[c] = cells[c].empty() ? 0 : 1;
    }
    glTextureSubImage3D(m_pathOccupancy, 0, 0, 0, 0, kPathR, kPathTheta, kPathPhi, GL_RED_INTEGER, GL_UNSIGNED_BYTE,
                        occupied.data());
    if (segments.empty()) segments.push_back({});  // buffers may not be empty
    if (indices.empty()) indices.push_back(0);
    glNamedBufferData(m_pathSegments, sizeof(GpuPathSegment) * segments.size(), segments.data(), GL_DYNAMIC_DRAW);
    glNamedBufferData(m_pathCells, sizeof(glm::uvec2) * ranges.size(), ranges.data(), GL_DYNAMIC_DRAW);
    glNamedBufferData(m_pathIndex, sizeof(uint32_t) * indices.size(), indices.data(), GL_DYNAMIC_DRAW);
    m_pathSegmentCount = static_cast<int>(indices.size());
}

void Renderer::shutdown() {
    m_hdr.destroy();
    m_ldr.destroy();
    for (auto& b : m_bloom) b.destroy();
    if (m_vao) glDeleteVertexArrays(1, &m_vao);
    if (m_bodyUbo) glDeleteBuffers(1, &m_bodyUbo);
    GLuint pathBuffers[] = {m_pathSegments, m_pathCells, m_pathIndex};
    glDeleteBuffers(3, pathBuffers);
    m_pathSegments = m_pathCells = m_pathIndex = 0;
    glDeleteTextures(1, &m_pathOccupancy);
    m_pathOccupancy = 0;
    if (m_outflowTexture) glDeleteTextures(1, &m_outflowTexture);
    m_outflowTexture = 0;
    m_outflowRows = 0;
    GLuint textures[] = {m_blackbodyTexture, m_skyTexture, m_placeholderTexture, m_placeholderVolume,
                         m_placeholderOccupancy, m_historyTexture};
    glDeleteTextures(6, textures);
    glDeleteQueries(kQueries, m_queries);
    m_vao = m_bodyUbo = m_blackbodyTexture = m_skyTexture = m_placeholderTexture = 0;
    m_placeholderVolume = m_placeholderOccupancy = m_historyTexture = 0;
}

void Renderer::uploadHistory(History& history) {
    constexpr int N = History::kColumns;
    constexpr int rows = History::kRows * History::kSlots;
    glPixelStorei(GL_UNPACK_ROW_LENGTH, N);

    // A new time step (the scene scale changed): everything was laid down again.
    if (history.generation != m_historyGeneration) {
        glTextureSubImage2D(m_historyTexture, 0, 0, 0, N, rows, GL_RGBA, GL_FLOAT, history.image.data());
        history.dirtyRows.clear();
        m_historyUploaded = history.newest;
        m_historyGeneration = history.generation;
        glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
        return;
    }

    // Whole rows for bodies whose past was just filled in.
    for (int slot : history.dirtyRows) {
        glTextureSubImage2D(m_historyTexture, 0, 0, History::kRows * slot, N, History::kRows, GL_RGBA, GL_FLOAT,
                            &history.image[static_cast<size_t>(History::kRows * slot) * N]);
    }
    history.dirtyRows.clear();

    // Columns recorded since the last upload, for every row (the ring may wrap).
    long long from = std::max(m_historyUploaded + 1, history.newest - N + 1);
    while (from <= history.newest) {
        int col = static_cast<int>(((from % N) + N) % N);
        int count = static_cast<int>(std::min<long long>(history.newest - from + 1, N - col));
        glTextureSubImage2D(m_historyTexture, 0, col, 0, count, rows, GL_RGBA, GL_FLOAT, &history.image[col]);
        from += count;
    }
    m_historyUploaded = history.newest;
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
}

void Renderer::drawFullscreen() const {
    glBindVertexArray(m_vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
}

void Renderer::readTimer() {
    // Read the query issued kQueries - 1 frames ago so we never stall the GPU.
    GLuint q = m_queries[(m_queryFrame + 1) % kQueries];
    if (m_queryFrame + 1 < kQueries) return;
    GLint available = 0;
    glGetQueryObjectiv(q, GL_QUERY_RESULT_AVAILABLE, &available);
    if (available) {
        GLuint64 ns = 0;
        glGetQueryObjectui64v(q, GL_QUERY_RESULT, &ns);
        float ms = static_cast<float>(ns * 1e-6);
        m_gpuMs = m_gpuMs == 0.0f ? ms : 0.9f * m_gpuMs + 0.1f * ms;
    }
}

void Renderer::render(const FrameInputs& in) {
    const Settings& s = *in.settings;
    const Camera& cam = *in.camera;

    const float scale = in.renderScale > 0.0f ? in.renderScale : s.renderScale;
    const int tw = std::max(1, static_cast<int>(std::lround(in.width * scale)));
    const int th = std::max(1, static_cast<int>(std::lround(in.height * scale)));
    m_hdr.ensure(tw, th, GL_RGBA16F);
    m_ldr.ensure(in.width, in.height, GL_RGBA8);
    m_bloomCount = 0;
    for (int i = 0, w = tw / 2, h = th / 2; i < kBloomLevels && w >= 2 && h >= 2; ++i, w /= 2, h /= 2) {
        m_bloom[i].ensure(w, h, GL_R11F_G11F_B10F);
        m_bloomCount = i + 1;
    }

    readTimer();
    glBeginQuery(GL_TIME_ELAPSED, m_queries[m_queryFrame % kQueries]);

    // --- 1. Ray trace ---------------------------------------------------------
    // Render units: lengths and times over the scene scale, mass 1 / scale.
    const double spin = in.spin;
    const double unit = in.scale;  // the render unit, in M
    const double inv = 1.0 / unit;
    const bool prograde = s.diskPrograde;
    const double rh = kerr::horizonRadius(spin);
    const CentralView fallback;
    const CentralView& cv = in.central ? *in.central : fallback;
    double diskInner = s.diskInnerAtIsco ? kerr::iscoRadius(spin, prograde) : std::max<double>(s.diskInner, rh + 0.05);
    double diskOuter = std::max(s.diskOuter, static_cast<float>(diskInner) + 0.5f);
    if (in.central) {
        diskInner = cv.diskInner;
        diskOuter = std::max(cv.diskOuter, cv.diskInner * 1.05);
    }
    const float tanHalfFov = std::tan(0.5f * glm::radians(cam.fovDeg));
    const float pixelAngle = 2.0f * tanHalfFov / static_cast<float>(th);

    // Compile in only what the scene needs (see traceProgram).
    bool lenses = false, planets = false;
    for (int i = 0; i < in.bodies->info.x; ++i) {
        const auto kind = static_cast<BodyKind>(std::lround(in.bodies->bodies[i].params.z));
        lenses = lenses || kind == BodyKind::MiniBlackHole;
        planets = planets || kind == BodyKind::RockyPlanet || kind == BodyKind::GasGiant || kind == BodyKind::Moon;
    }
    // The Sun's 30.4 nm map: all of it the first time, then a band of rows a frame.
    if (cv.sdoMap) {
        bakeSdo(cv, !m_sdoReady);
        m_sdoReady = true;
    } else {
        m_sdoReady = false;
    }
    const FlareView noFlare;
    const FlareView& fv = in.flare ? *in.flare : noFlare;
    const bool outflow = fv.active && cv.surface <= 0.0 && !fv.table.empty();
    ShaderProgram& trace = traceProgram(cv.surface > 0.0 ? std::max(cv.look, 1) : 0, lenses, planets, outflow);
    glBindFramebuffer(GL_FRAMEBUFFER, m_hdr.fbo);
    glViewport(0, 0, tw, th);
    trace.use();
    trace.set("uResolution", glm::vec2(tw, th));
    trace.set("uSpin", static_cast<float>(spin * inv));
    trace.set("uMass", static_cast<float>(inv));
    trace.set("uCamPos", glm::vec3(static_cast<float>(cam.r * inv), static_cast<float>(cam.theta),
                                     static_cast<float>(cam.phi)));
    trace.set("uCamBasis", cam.basis());
    {
        // Each frame vector as (u_t, u_phi, u_theta, u^r); see uCamPhoton in trace.frag.
        // Computed in M (double precision), then u_phi and u_theta, which carry a length,
        // go to render units; u_t and u^r are dimensionless.
        const kerr::Frame f = cam.frame(spin);
        const kerr::MetricTensor g = kerr::metricTensor(spin, cam.r, cam.theta);
        glm::mat4 m(0.0f);
        for (int j = 0; j < 4; ++j) {
            const kerr::Vec4 low = kerr::lower(g, f.e[j]);
            m[j] = glm::vec4(static_cast<float>(low[0]), static_cast<float>(low[3] * inv),
                             static_cast<float>(low[2] * inv), static_cast<float>(f.e[j][1]));
        }
        trace.set("uCamPhoton", m);
    }
    trace.set("uTanHalfFov", tanHalfFov);
    trace.set("uPixelAngle", pixelAngle);
    trace.set("uMaxSteps", s.maxSteps);
    trace.set("uStepSize", s.stepSize);
    // Rays must get past everything that glows: the bodies, the far debris, the outflow.
    float escape = std::max({s.escapeRadius, static_cast<float>(cam.r * inv) * 1.5f,
                             static_cast<float>(in.reach * inv) * 1.2f});
    if (in.particles && in.particles->active()) escape = std::max(escape, in.particles->gridRMax() * 1.02f);
    if (outflow) escape = std::max(escape, static_cast<float>(fv.extent * inv) * 1.05f + 5.0f);
    trace.set("uEscapeRadius", escape);

    // With the disk off (a quiet hole), a tidal disruption's returning gas forms its own
    // small disk, which shines as the flare does.
    const bool flareDisk = !s.diskEnabled && fv.diskGlow > 1e-4f && cv.surface <= 0.0;
    if (flareDisk) {
        diskInner = kerr::iscoRadius(spin, true);
        diskOuter = std::max(fv.diskOuter, diskInner * 1.5);
    }
    trace.set("uDiskEnabled", s.diskEnabled || flareDisk ? 1 : 0);
    trace.set("uDiskInner", static_cast<float>(diskInner * inv));
    trace.set("uDiskOuter", static_cast<float>(diskOuter * inv));
    trace.set("uDiskTemp", flareDisk ? fv.diskTemperature : s.diskTemperature);
    trace.set("uDiskBrightness", flareDisk ? s.diskBrightness * fv.diskGlow : s.diskBrightness);
    trace.set("uDiskOpacity", flareDisk ? s.diskOpacity * std::min(1.0f, 4.0f * fv.diskGlow) : s.diskOpacity);
    trace.set("uDiskTurbulence", s.diskTurbulence);
    trace.set("uDiskDirection", prograde || flareDisk ? 1.0f : -1.0f);
    // The disk's pattern is renewed every few dynamical times: kFlowPeriod M around a
    // hole, sqrt(scale) times that in render units around a bigger object.
    const double flowPeriod = kFlowPeriod * std::sqrt(unit);
    trace.set("uFlowTime", static_cast<float>(std::fmod(in.time * inv, flowPeriod)));
    trace.set("uFlowPeriod", static_cast<float>(flowPeriod));
    trace.set("uRelativity", s.relativity);

    const bool debris = in.particles && in.particles->active();
    trace.set("uDebrisEnabled", debris ? 1 : 0);
    trace.set("uGridLogR", in.particles ? in.particles->gridLogR()
                                        : glm::vec2(std::log(ParticleSystem::kGridRMin),
                                                    std::log(ParticleSystem::kGridRMax)));
    trace.set("uDebrisBrightness", s.debrisBrightness);
    trace.set("uDebrisOpacity", s.debrisOpacity);
    trace.set("uDebrisFilaments", s.debrisFilaments ? 1.0f : 0.0f);
    trace.set("uDebrisGlowFloor", s.debrisGlowFloor);

    trace.set("uPathsEnabled", m_pathSegmentCount > 0 ? 1 : 0);
    trace.set("uPathGridDims", glm::ivec3(kPathR, kPathTheta, kPathPhi));
    trace.set("uPathLogR", glm::vec2(std::log(ParticleSystem::kGridRMin), std::log(m_pathRMax)));
    trace.set("uPathWidth", s.pathWidth);
    trace.set("uSkipBody", in.skipBody);
    trace.set("uHistoryScale", static_cast<float>(inv));
    trace.set("uJetsEnabled", s.jetsEnabled ? 1 : 0);
    trace.set("uJetPower", s.jetPower);
    trace.set("uJetSpeed", std::clamp(s.jetSpeed, 0.0f, 0.995f));
    trace.set("uJetOpening", std::tan(glm::radians(s.jetOpeningDeg)));
    trace.set("uTime", static_cast<float>(std::fmod(in.time * inv, 10000.0)));

    // The central object.
    trace.set("uCentralKind", cv.kind);
    trace.set("uSurfaceR", static_cast<float>(cv.surface * inv));
    trace.set("uSurfaceOmega", static_cast<float>(cv.omega * unit));
    trace.set("uSurfacePhase", static_cast<float>(std::fmod(cv.phase, 2.0 * kerr::kPi)));
    trace.set("uSurfaceTurns", static_cast<float>(cv.turns));
    trace.set("uSurfaceShift", static_cast<float>(cv.azimuthShift));
    trace.set("uSurfaceTemp", cv.temperature);
    trace.set("uHotSpotTemp", cv.hotSpotTemperature);
    trace.set("uSurfaceBrightness", cv.brightness);
    trace.set("uSurfaceDetail", cv.detail ? 1 : 0);
    trace.set("uMagneticAxis", cv.magneticAxis);
    trace.set("uPulsar", cv.pulsar ? 1 : 0);
    trace.set("uMagnetar", cv.magnetar ? 1 : 0);
    trace.set("uSurfaceLook", cv.look);
    trace.set("uSurfaceKm", static_cast<float>(std::max(cv.radiusKm, 1e-3)));
    trace.set("uPlanetSeed", cv.planetSeed);
    trace.set("uPlanetLight", cv.planetLight);
    trace.set("uAtmosphere", cv.atmosphere ? 1 : 0);
    trace.set("uCoronaScale", cv.coronaScale);
    trace.set("uChromosphereScale", cv.chromosphereScale);
    trace.set("uChromosphereHeight", cv.chromosphereHeight);
    trace.set("uFilmColor", cv.filmColor ? 1 : 0);
    trace.set("uEnvelope", cv.envelope);
    trace.set("uFlowGlow", cv.flowGlow);
    trace.set("uGranuleFreq", cv.granuleFrequency);
    trace.set("uGranuleContrast", cv.granuleContrast);
    trace.set("uGranuleLane", cv.granuleLane);
    trace.set("uGranulePhase", cv.granulePhase);
    trace.set("uGranuleRate", static_cast<float>(cv.granuleRate * unit));  // per unit render time
    trace.set("uCapCos", cv.capCos);
    {
        glm::vec4 spots[32] = {};
        float umbra[32] = {};
        const int count = static_cast<int>(std::min<size_t>(cv.spots.size(), 32));
        for (int k = 0; k < count; ++k) {
            spots[k] = cv.spots[k];
            umbra[k] = k < static_cast<int>(cv.spotUmbra.size()) ? cv.spotUmbra[k] : 0.4f;
        }
        trace.set("uSpotCount", count);
        glProgramUniform4fv(trace.id(), trace.location("uSpots"), 32, &spots[0].x);
        glProgramUniform1fv(trace.id(), trace.location("uSpotUmbra"), 32, umbra);
        trace.set("uSpotDrop", cv.spotDrop);
        glm::vec4 regions[8] = {};
        const int regionCount = static_cast<int>(std::min<size_t>(cv.regions.size(), 8));
        for (int k = 0; k < regionCount; ++k) regions[k] = cv.regions[k];
        trace.set("uRegionCount", regionCount);
        glProgramUniform4fv(trace.id(), trace.location("uRegions"), 8, &regions[0].x);
        glm::vec4 flares[4] = {};
        float flareAges[4] = {};
        const int flareCount = static_cast<int>(std::min<size_t>(cv.flares.size(), 4));
        for (int k = 0; k < flareCount; ++k) {
            flares[k] = cv.flares[k];
            flareAges[k] = static_cast<float>(cv.flareAges[k] * inv);
        }
        trace.set("uFlareCount", flareCount);
        glProgramUniform4fv(trace.id(), trace.location("uFlares"), 4, &flares[0].x);
        glProgramUniform1fv(trace.id(), trace.location("uFlareAges"), 4, flareAges);
        trace.set("uFlareRise", static_cast<float>(std::max(cv.flareRise * inv, 1e-6)));
        trace.set("uFlareDecay", static_cast<float>(std::max(cv.flareDecay * inv, 1e-6)));
    }
    trace.set("uFunnelRadius", static_cast<float>(cv.funnelRadius * inv));
    trace.set("uFieldGlow", cv.fieldGlow);
    trace.set("uFieldTwist", cv.fieldTwist);
    trace.set("uFieldReach", static_cast<float>(cv.fieldReach * inv));
    trace.set("uPairGlow", cv.pairGlow);
    trace.set("uBeamGlow", cv.beamGlow);
    trace.set("uBeacon", cv.beacon ? 1 : 0);
    trace.set("uHaloGlow", cv.haloGlow);
    {
        glm::vec4 quakes[4] = {}, axes[4] = {};
        const int count = static_cast<int>(std::min<size_t>(std::min(cv.quakes.size(), cv.quakeAxes.size()), 4));
        for (int k = 0; k < count; ++k) {
            // Into render units: ages and decay times over the scale, frequency times it.
            quakes[k] = glm::vec4(static_cast<float>(cv.quakes[k].x * inv), cv.quakes[k].y,
                                  static_cast<float>(cv.quakes[k].z * unit), static_cast<float>(cv.quakes[k].w * inv));
            axes[k] = cv.quakeAxes[k];
        }
        trace.set("uQuakeCount", count);
        glProgramUniform4fv(trace.id(), trace.location("uQuakes"), 4, &quakes[0].x);
        glProgramUniform4fv(trace.id(), trace.location("uQuakeAxes"), 4, &axes[0].x);
    }
    trace.set("uBeamWidth", std::max(cv.beamWidth, 1e-3f));
    trace.set("uCompanion", glm::vec4(glm::vec3(cv.companion) * static_cast<float>(inv),
                                      cv.companion.w * static_cast<float>(inv)));
    trace.set("uCompanionAxis", cv.companionAxis);
    trace.set("uCompanionTemp", cv.companionTemperature);
    {
        glm::vec3 stream[16] = {};
        const int count = static_cast<int>(std::min<size_t>(cv.stream.size(), 16));
        for (int k = 0; k < count; ++k) stream[k] = cv.stream[k] * static_cast<float>(inv);
        trace.set("uStreamCount", count);
        glProgramUniform3fv(trace.id(), trace.location("uStream"), 16, &stream[0].x);
    }
    trace.set("uStreamWidth", static_cast<float>(cv.streamWidth * inv));
    trace.set("uBinaryPhase", cv.binaryPhase);
    trace.set("uSpotPhase", cv.spotPhase);
    trace.set("uSunDirection", glm::normalize(cv.sunDirection));
    trace.set("uBodyLight", cv.bodyLight);
    trace.set("uBodyLightTemp", cv.bodyLightTemperature);
    // Light falling off as 1/d^2: the scale carries a length squared.
    trace.set("uBodyLightScale", static_cast<float>(cv.bodyLight == 1 ? cv.bodyLightScale * inv * inv
                                                                       : cv.bodyLightScale));
    {
        const int count = static_cast<int>(std::min<size_t>(cv.impacts.size(), 8));
        glm::vec4 impacts[8] = {};
        float ages[8] = {};
        for (int k = 0; k < count; ++k) {
            impacts[k] = cv.impacts[k];
            ages[k] = static_cast<float>(std::min(cv.impactAges[k] * inv, 1e30));
        }
        trace.set("uImpactCount", count);
        glProgramUniform4fv(trace.id(), trace.location("uImpacts"), 8, &impacts[0].x);
        glProgramUniform1fv(trace.id(), trace.location("uImpactAges"), 8, ages);
        trace.set("uImpactFade", static_cast<float>(cv.impactFade * inv));
        trace.set("uFlash", glm::vec3(static_cast<float>(std::min(cv.flashAge * inv, 1e30)),
                                        static_cast<float>(cv.flashStrength),
                                        static_cast<float>(std::max(cv.flashFade * inv, 1e-6))));
    }
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, m_pathSegments);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, m_pathCells);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, m_pathIndex);

    if (outflow) {
        // The history table (render units: times over the scale).
        const int rows = static_cast<int>(fv.table.size());
        if (rows != m_outflowRows) {
            if (m_outflowTexture) glDeleteTextures(1, &m_outflowTexture);
            glCreateTextures(GL_TEXTURE_1D, 1, &m_outflowTexture);
            glTextureStorage1D(m_outflowTexture, 1, GL_RGBA32F, rows);
            glTextureParameteri(m_outflowTexture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTextureParameteri(m_outflowTexture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTextureParameteri(m_outflowTexture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            m_outflowRows = rows;
        }
        glTextureSubImage1D(m_outflowTexture, 0, 0, rows, GL_RGBA, GL_FLOAT, fv.table.data());
        glm::vec3 axis = glm::normalize(fv.axis);
        glm::vec3 e1 = glm::normalize(glm::cross(axis, std::abs(axis.z) < 0.9f ? glm::vec3(0, 0, 1) : glm::vec3(1, 0, 0)));
        glm::vec3 e2 = glm::cross(axis, e1);
        trace.set("uOutflowAxis", axis);
        trace.set("uOutflowE1", e1);
        trace.set("uOutflowE2", e2);
        trace.set("uOutflowClock", static_cast<float>(fv.clock * inv));
        trace.set("uOutflowTable", glm::vec2(static_cast<float>(fv.tableStart * inv),
                                             static_cast<float>(unit / std::max(fv.tableStep, 1e-9))));
        trace.set("uOutflowRows", static_cast<float>(rows));
        trace.set("uOutflowStart", static_cast<float>(fv.start * inv));
        trace.set("uOutflowSpeed", static_cast<float>(std::clamp(fv.speed, 0.005, 0.9)));
        trace.set("uOutflowLaunch", static_cast<float>(fv.launchRadius * inv));
        trace.set("uOutflowExtent", static_cast<float>(fv.extent * inv));
        trace.set("uOutflowTanMax", std::tan(static_cast<float>(std::min(fv.maxOpening * 1.6, 1.45))));
        trace.set("uOutflowMaxAngle", static_cast<float>(std::min(fv.maxOpening * 1.6, 1.5)));
        {
            float cone[FlareView::kConeRows] = {};
            for (int k = 0; k < FlareView::kConeRows && k < static_cast<int>(fv.cone.size()); ++k) cone[k] = fv.cone[k];
            glProgramUniform1fv(trace.id(), trace.location("uOutflowCone"), FlareView::kConeRows, cone);
            const double logR0 = std::log(fv.launchRadius * inv);
            const double span = std::max(std::log(fv.extent / fv.launchRadius), 1e-3);
            trace.set("uOutflowConeR", glm::vec2(static_cast<float>(logR0),
                                                 static_cast<float>((FlareView::kConeRows - 1) / span)));
        }
        trace.set("uOutflowThickness", fv.thickness);
        trace.set("uOutflowLight", fv.brightness);
        trace.set("uOutflowFalloff", fv.falloff);
        trace.set("uFlareTemp", fv.flareTemperature);
        trace.set("uOutflowGlow", fv.glow);
        trace.set("uOutflowEnvelope", fv.envelope);
        trace.set("uOutflowSwirl", fv.swirl);
        glBindTextureUnit(12, m_outflowTexture);
    }

    const int skyMode = (s.skyMode == 2 && !hasSkyImage()) ? 0 : s.skyMode;
    trace.set("uSkyMode", skyMode);
    trace.set("uStarBrightness", s.starBrightness);
    trace.set("uGalaxyBrightness", s.galaxyBrightness);
    trace.set("uSkyImageBrightness", s.skyImageBrightness);
    trace.set("uSkyLod", m_skyWidth > 0
                               ? std::max(0.0f, std::log2(pixelAngle * m_skyWidth / 6.2831853f))
                               : 0.0f);
    trace.set("uBlackbodyRange", glm::vec2(static_cast<float>(blackbody::kLogTMin),
                                             static_cast<float>(blackbody::kLogTMax)));
    trace.set("uDebugView", s.debugView);

    glNamedBufferSubData(m_bodyUbo, 0, sizeof(GpuBodies), in.bodies);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, m_bodyUbo);
    glBindTextureUnit(0, m_blackbodyTexture);
    glBindTextureUnit(1, m_skyTexture ? m_skyTexture : m_placeholderTexture);
    glBindTextureUnit(2, m_historyTexture);
    glBindTextureUnit(3, debris ? in.particles->densityTexture() : m_placeholderVolume);
    glBindTextureUnit(4, debris ? in.particles->velocityTexture() : m_placeholderVolume);
    glBindTextureUnit(5, debris ? in.particles->occupancyTexture() : m_placeholderOccupancy);
    glBindTextureUnit(6, m_pathOccupancy);
    const GLuint impactMap = in.particles ? in.particles->impactTexture() : 0;
    glBindTextureUnit(7, impactMap ? impactMap : m_placeholderTexture);
    glBindTextureUnit(11, cv.sdoMap && m_sdoMap ? m_sdoMap : m_placeholderTexture);
    {
        const bool earth = cv.look == 5 && cv.planetMaps && m_earthDay && m_earthClouds && m_earthNight;
        const bool jupiter = cv.look == 6 && cv.planetMaps && m_jupiterMap;
        glBindTextureUnit(8, earth ? m_earthDay : (jupiter ? m_jupiterMap : m_placeholderTexture));
        glBindTextureUnit(9, earth ? m_earthClouds : m_placeholderTexture);
        glBindTextureUnit(10, earth ? m_earthNight : m_placeholderTexture);
        trace.set("uPlanetMaps", earth || jupiter ? 1 : 0);
    }
    drawFullscreen();

    // --- 2. Bloom: progressive downsample, then tent-filtered upsample --------
    // Debug views are shown raw: no bloom, exposure or tone curve.
    const bool debugView = s.debugView != 0;
    const bool bloom = s.bloomStrength > 0.0f && m_bloomCount > 0 && !debugView;
    if (bloom) {
        m_bloomDown.use();
        GLuint src = m_hdr.texture;
        glm::vec2 srcSize(tw, th);
        for (int i = 0; i < m_bloomCount; ++i) {
            glBindFramebuffer(GL_FRAMEBUFFER, m_bloom[i].fbo);
            glViewport(0, 0, m_bloom[i].width, m_bloom[i].height);
            m_bloomDown.set("uSrcTexel", 1.0f / srcSize);
            m_bloomDown.set("uFirst", i == 0 ? 1 : 0);
            glBindTextureUnit(0, src);
            drawFullscreen();
            src = m_bloom[i].texture;
            srcSize = glm::vec2(m_bloom[i].width, m_bloom[i].height);
        }

        m_bloomUp.use();
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
        for (int i = m_bloomCount - 1; i > 0; --i) {
            glBindFramebuffer(GL_FRAMEBUFFER, m_bloom[i - 1].fbo);
            glViewport(0, 0, m_bloom[i - 1].width, m_bloom[i - 1].height);
            m_bloomUp.set("uSrcTexel", glm::vec2(1.0f / m_bloom[i].width, 1.0f / m_bloom[i].height));
            glBindTextureUnit(0, m_bloom[i].texture);
            drawFullscreen();
        }
        glDisable(GL_BLEND);
    }

    // --- 3. Tone map to 8-bit ---------------------------------------------------
    glBindFramebuffer(GL_FRAMEBUFFER, m_ldr.fbo);
    glViewport(0, 0, in.width, in.height);
    m_tonemap.use();
    m_tonemap.set("uExposure", debugView ? 1.0f : s.exposure * in.exposureBoost);
    m_tonemap.set("uBloomStrength", bloom ? s.bloomStrength : 0.0f);
    m_tonemap.set("uBloomLevels", static_cast<float>(std::max(m_bloomCount, 1)));
    m_tonemap.set("uTonemapper", debugView ? 2 : s.tonemapper);
    glBindTextureUnit(0, m_hdr.texture);
    glBindTextureUnit(1, bloom ? m_bloom[0].texture : m_hdr.texture);
    drawFullscreen();

    glEndQuery(GL_TIME_ELAPSED);
    ++m_queryFrame;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Renderer::present(int width, int height) const {
    // Fit the image inside the window, keeping its aspect ratio.
    const float imageAspect = static_cast<float>(m_ldr.width) / static_cast<float>(m_ldr.height);
    const float windowAspect = static_cast<float>(width) / static_cast<float>(height);
    int w = width, h = height, x = 0, y = 0;
    if (imageAspect > windowAspect * 1.001f) {
        h = static_cast<int>(width / imageAspect);
        y = (height - h) / 2;
    } else if (imageAspect < windowAspect * 0.999f) {
        w = static_cast<int>(height * imageAspect);
        x = (width - w) / 2;
    }
    if (w != width || h != height) {
        const float black[4] = {0.0f, 0.0f, 0.0f, 1.0f};
        glClearNamedFramebufferfv(0, GL_COLOR, 0, black);
    }
    glBlitNamedFramebuffer(m_ldr.fbo, 0, 0, 0, m_ldr.width, m_ldr.height, x, y, x + w, y + h,
                           GL_COLOR_BUFFER_BIT, GL_LINEAR);
}

void Renderer::readOutput(std::vector<unsigned char>& bgra, int& width, int& height) const {
    width = m_ldr.width;
    height = m_ldr.height;
    const size_t rowBytes = static_cast<size_t>(width) * 4;
    std::vector<unsigned char> bottomUp(rowBytes * height);
    glGetTextureImage(m_ldr.texture, 0, GL_BGRA, GL_UNSIGNED_BYTE, static_cast<GLsizei>(bottomUp.size()),
                      bottomUp.data());
    bgra.resize(bottomUp.size());
    for (int row = 0; row < height; ++row) {
        std::copy_n(&bottomUp[rowBytes * (height - 1 - row)], rowBytes, &bgra[rowBytes * row]);
    }
}

bool Renderer::saveScreenshot(const std::string& path, std::string& error) const {
    if (!m_ldr.texture) {
        error = "nothing rendered yet";
        return false;
    }
    std::vector<unsigned char> pixels(static_cast<size_t>(m_ldr.width) * m_ldr.height * 4);
    glGetTextureImage(m_ldr.texture, 0, GL_RGBA, GL_UNSIGNED_BYTE, static_cast<GLsizei>(pixels.size()),
                      pixels.data());
    for (size_t i = 3; i < pixels.size(); i += 4) pixels[i] = 255;
    stbi_flip_vertically_on_write(1);
    if (!stbi_write_png(path.c_str(), m_ldr.width, m_ldr.height, 4, pixels.data(), m_ldr.width * 4)) {
        error = "could not write " + path;
        return false;
    }
    return true;
}

namespace {
// An equirectangular map (row 0 north after the flip below: v = 0 is the south pole),
// mipmapped and wrapping in longitude; 0 if the file is missing.
GLuint loadMap(const std::filesystem::path& path, bool srgb) {
    stbi_set_flip_vertically_on_load(1);
    int w = 0, h = 0, channels = 0;
    unsigned char* data = stbi_load(path.string().c_str(), &w, &h, &channels, 4);
    if (!data) return 0;
    GLuint texture = 0;
    glCreateTextures(GL_TEXTURE_2D, 1, &texture);
    glTextureStorage2D(texture, 1 + static_cast<int>(std::log2(std::max(w, h))), srgb ? GL_SRGB8_ALPHA8 : GL_RGBA8, w,
                       h);
    glTextureSubImage2D(texture, 0, 0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, data);
    stbi_image_free(data);
    glGenerateTextureMipmap(texture);
    glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    return texture;
}
}  // namespace

void Renderer::bakeSdo(const CentralView& cv, bool all) {
    if (!m_sdoBake.valid()) return;
    if (!m_sdoMap) {
        glCreateTextures(GL_TEXTURE_2D, 1, &m_sdoMap);
        glTextureStorage2D(m_sdoMap, 1 + static_cast<int>(std::log2(kSdoWidth)), GL_R16F, kSdoWidth, kSdoHeight);
        glTextureParameteri(m_sdoMap, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTextureParameteri(m_sdoMap, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTextureParameteri(m_sdoMap, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTextureParameteri(m_sdoMap, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glCreateFramebuffers(1, &m_sdoFbo);
        glNamedFramebufferTexture(m_sdoFbo, GL_COLOR_ATTACHMENT0, m_sdoMap, 0);
    }
    m_sdoBake.use();
    glm::vec4 regions[8] = {};
    const int count = static_cast<int>(std::min<size_t>(cv.regions.size(), 8));
    for (int k = 0; k < count; ++k) regions[k] = cv.regions[k];
    m_sdoBake.set("uRegionCount", count);
    glProgramUniform4fv(m_sdoBake.id(), m_sdoBake.location("uRegions"), 8, &regions[0].x);
    m_sdoBake.set("uGranuleFreq", cv.granuleFrequency);
    m_sdoBake.set("uTime", cv.granulePhase * 0.03f);
    m_sdoBake.set("uMapSize", glm::vec2(kSdoWidth, kSdoHeight));
    glBindFramebuffer(GL_FRAMEBUFFER, m_sdoFbo);
    glViewport(0, 0, kSdoWidth, kSdoHeight);
    glEnable(GL_SCISSOR_TEST);
    const int rows = kSdoHeight / kSdoBands;
    glScissor(0, all ? 0 : m_sdoBand * rows, kSdoWidth, all ? kSdoHeight : rows);
    glBindVertexArray(m_vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glDisable(GL_SCISSOR_TEST);
    glGenerateTextureMipmap(m_sdoMap);
    m_sdoBand = (m_sdoBand + 1) % kSdoBands;
}

void Renderer::loadPlanetMaps() {
    // Next to the shaders' folder: the source tree while developing, the executable's
    // folder otherwise (both hold assets/planets).
    const std::filesystem::path dir = m_shaderDir.parent_path() / "assets" / "planets";
    m_earthDay = loadMap(dir / "earth_day.jpg", true);
    m_earthClouds = loadMap(dir / "earth_clouds.jpg", false);
    m_earthNight = loadMap(dir / "earth_night.jpg", false);
    m_jupiterMap = loadMap(dir / "jupiter.jpg", true);
    std::printf("Planet maps: %s\n", m_earthDay && m_earthClouds && m_earthNight && m_jupiterMap
                                         ? "loaded"
                                         : ("not found in " + dir.string()).c_str());
}

bool Renderer::loadSkyImage(const std::string& path, std::string& error) {
    stbi_set_flip_vertically_on_load(1);
    int w = 0, h = 0, channels = 0;
    GLuint texture = 0;
    if (stbi_is_hdr(path.c_str())) {
        float* data = stbi_loadf(path.c_str(), &w, &h, &channels, 3);
        if (!data) {
            error = stbi_failure_reason();
            return false;
        }
        glCreateTextures(GL_TEXTURE_2D, 1, &texture);
        glTextureStorage2D(texture, 1 + static_cast<int>(std::log2(std::max(w, h))), GL_RGB16F, w, h);
        glTextureSubImage2D(texture, 0, 0, 0, w, h, GL_RGB, GL_FLOAT, data);
        stbi_image_free(data);
    } else {
        unsigned char* data = stbi_load(path.c_str(), &w, &h, &channels, 4);
        if (!data) {
            error = stbi_failure_reason();
            return false;
        }
        glCreateTextures(GL_TEXTURE_2D, 1, &texture);
        glTextureStorage2D(texture, 1 + static_cast<int>(std::log2(std::max(w, h))), GL_SRGB8_ALPHA8, w, h);
        glTextureSubImage2D(texture, 0, 0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, data);
        stbi_image_free(data);
    }
    glGenerateTextureMipmap(texture);
    glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    if (m_skyTexture) glDeleteTextures(1, &m_skyTexture);
    m_skyTexture = texture;
    m_skyWidth = w;
    return true;
}

void Renderer::hotReload() {
    for (ShaderProgram& p : m_trace) p.reloadIfChanged();
    m_bloomDown.reloadIfChanged();
    m_sdoBake.reloadIfChanged();
    m_bloomUp.reloadIfChanged();
    m_tonemap.reloadIfChanged();
}

std::string Renderer::shaderErrors() const {
    std::string all;
    for (const ShaderProgram& p : m_trace) {
        if (!p.error().empty()) all += p.error();
    }
    for (const ShaderProgram* p : {&m_bloomDown, &m_bloomUp, &m_tonemap}) {
        if (!p->error().empty()) all += p->error();
    }
    return all;
}
