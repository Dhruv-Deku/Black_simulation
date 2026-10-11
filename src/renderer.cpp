#include "renderer.h"

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

// Cells binning path segments for the ray tracer, over the debris grid's r range.
constexpr int kPathR = 64, kPathTheta = 48, kPathPhi = 96;
constexpr int kPathCells = kPathR * kPathTheta * kPathPhi;

struct GpuPathSegment {
    glm::vec4 a, b, color, dash;
};
static_assert(sizeof(GpuPathSegment) == 64, "std430 layout");

// The path cell containing a Cartesian point, or -1 outside the grid.
int pathCell(const glm::vec3& p, double spin) {
    auto [r, theta, phi] = kerr::fromCartesian(spin, p.x, p.y, p.z);
    const double lr = (std::log(r) - std::log(ParticleSystem::kGridRMin)) /
                      (std::log(ParticleSystem::kGridRMax) - std::log(ParticleSystem::kGridRMin));
    if (lr < 0.0 || lr >= 1.0) return -1;
    double f = phi / (2.0 * kerr::kPi);
    f -= std::floor(f);
    const int x = std::min(static_cast<int>(lr * kPathR), kPathR - 1);
    const int y = std::clamp(static_cast<int>(theta / kerr::kPi * kPathTheta), 0, kPathTheta - 1);
    const int z = std::min(static_cast<int>(f * kPathPhi), kPathPhi - 1);
    return (z * kPathTheta + y) * kPathR + x;
}
}  // namespace

bool Renderer::init(const std::filesystem::path& shaderDir, std::string& error) {
    bool ok = m_trace.build(shaderDir, "fullscreen.vert", "trace.frag");
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
    setPaths({}, 0.0);
    return ok;
}

void Renderer::setPaths(const std::vector<PathLine>& lines, double spin) {
    std::vector<GpuPathSegment> segments;
    std::vector<std::vector<uint32_t>> cells(kPathCells);
    for (const PathLine& line : lines) {
        const size_t n = line.points.size();
        if (n < 2) continue;
        float along = 0.0f;
        int prevCell = pathCell(line.points[0], spin);
        for (size_t i = 1; i < n; ++i) {
            const glm::vec3& a = line.points[i - 1];
            const glm::vec3& b = line.points[i];
            const float len = glm::length(b - a);
            const float fa = static_cast<float>(i - 1) / static_cast<float>(n - 1);
            const float fb = static_cast<float>(i) / static_cast<float>(n - 1);
            const int cell = pathCell(b, spin);
            const int mid = pathCell(0.5f * (a + b), spin);
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
    const double spin = in.spin;
    const bool prograde = s.diskPrograde;
    const double rh = kerr::horizonRadius(spin);
    const double diskInner = s.diskInnerAtIsco ? kerr::iscoRadius(spin, prograde)
                                               : std::max<double>(s.diskInner, rh + 0.05);
    const float tanHalfFov = std::tan(0.5f * glm::radians(cam.fovDeg));
    const float pixelAngle = 2.0f * tanHalfFov / static_cast<float>(th);

    glBindFramebuffer(GL_FRAMEBUFFER, m_hdr.fbo);
    glViewport(0, 0, tw, th);
    m_trace.use();
    m_trace.set("uResolution", glm::vec2(tw, th));
    m_trace.set("uSpin", static_cast<float>(spin));
    m_trace.set("uCamPos", glm::vec3(static_cast<float>(cam.r), static_cast<float>(cam.theta),
                                     static_cast<float>(cam.phi)));
    m_trace.set("uCamBasis", cam.basis());
    {
        // Each frame vector as (u_t, u_phi, u_theta, u^r); see uCamPhoton in trace.frag.
        const kerr::Frame f = cam.frame(spin);
        const kerr::MetricTensor g = kerr::metricTensor(spin, cam.r, cam.theta);
        glm::mat4 m(0.0f);
        for (int j = 0; j < 4; ++j) {
            const kerr::Vec4 low = kerr::lower(g, f.e[j]);
            m[j] = glm::vec4(static_cast<float>(low[0]), static_cast<float>(low[3]), static_cast<float>(low[2]),
                             static_cast<float>(f.e[j][1]));
        }
        m_trace.set("uCamPhoton", m);
    }
    m_trace.set("uTanHalfFov", tanHalfFov);
    m_trace.set("uPixelAngle", pixelAngle);
    m_trace.set("uMaxSteps", s.maxSteps);
    m_trace.set("uStepSize", s.stepSize);
    m_trace.set("uEscapeRadius", std::max(s.escapeRadius, static_cast<float>(cam.r) * 1.5f));

    m_trace.set("uDiskEnabled", s.diskEnabled ? 1 : 0);
    m_trace.set("uDiskInner", static_cast<float>(diskInner));
    m_trace.set("uDiskOuter", std::max(s.diskOuter, static_cast<float>(diskInner) + 0.5f));
    m_trace.set("uDiskTemp", s.diskTemperature);
    m_trace.set("uDiskBrightness", s.diskBrightness);
    m_trace.set("uDiskOpacity", s.diskOpacity);
    m_trace.set("uDiskTurbulence", s.diskTurbulence);
    m_trace.set("uDiskDirection", prograde ? 1.0f : -1.0f);
    m_trace.set("uFlowTime", static_cast<float>(std::fmod(in.time, kFlowPeriod)));
    m_trace.set("uFlowPeriod", static_cast<float>(kFlowPeriod));
    m_trace.set("uRelativity", s.relativity);

    const bool debris = in.particles && in.particles->active();
    m_trace.set("uDebrisEnabled", debris ? 1 : 0);
    m_trace.set("uGridLogR", glm::vec2(std::log(ParticleSystem::kGridRMin), std::log(ParticleSystem::kGridRMax)));
    m_trace.set("uDebrisBrightness", s.debrisBrightness);
    m_trace.set("uDebrisOpacity", s.debrisOpacity);

    m_trace.set("uPathsEnabled", m_pathSegmentCount > 0 ? 1 : 0);
    m_trace.set("uPathGridDims", glm::ivec3(kPathR, kPathTheta, kPathPhi));
    m_trace.set("uPathWidth", s.pathWidth);
    m_trace.set("uSkipBody", in.skipBody);
    m_trace.set("uJetsEnabled", s.jetsEnabled ? 1 : 0);
    m_trace.set("uJetPower", s.jetPower);
    m_trace.set("uJetSpeed", std::clamp(s.jetSpeed, 0.0f, 0.995f));
    m_trace.set("uJetOpening", std::tan(glm::radians(s.jetOpeningDeg)));
    m_trace.set("uTime", static_cast<float>(std::fmod(in.time, 10000.0)));
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, m_pathSegments);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, m_pathCells);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, m_pathIndex);

    const int skyMode = (s.skyMode == 2 && !hasSkyImage()) ? 0 : s.skyMode;
    m_trace.set("uSkyMode", skyMode);
    m_trace.set("uStarBrightness", s.starBrightness);
    m_trace.set("uGalaxyBrightness", s.galaxyBrightness);
    m_trace.set("uSkyImageBrightness", s.skyImageBrightness);
    m_trace.set("uSkyLod", m_skyWidth > 0
                               ? std::max(0.0f, std::log2(pixelAngle * m_skyWidth / 6.2831853f))
                               : 0.0f);
    m_trace.set("uBlackbodyRange", glm::vec2(static_cast<float>(blackbody::kLogTMin),
                                             static_cast<float>(blackbody::kLogTMax)));
    m_trace.set("uDebugView", s.debugView);

    glNamedBufferSubData(m_bodyUbo, 0, sizeof(GpuBodies), in.bodies);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, m_bodyUbo);
    glBindTextureUnit(0, m_blackbodyTexture);
    glBindTextureUnit(1, m_skyTexture ? m_skyTexture : m_placeholderTexture);
    glBindTextureUnit(2, m_historyTexture);
    glBindTextureUnit(3, debris ? in.particles->densityTexture() : m_placeholderVolume);
    glBindTextureUnit(4, debris ? in.particles->velocityTexture() : m_placeholderVolume);
    glBindTextureUnit(5, debris ? in.particles->occupancyTexture() : m_placeholderOccupancy);
    glBindTextureUnit(6, m_pathOccupancy);
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
    m_tonemap.set("uExposure", debugView ? 1.0f : s.exposure);
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
    m_trace.reloadIfChanged();
    m_bloomDown.reloadIfChanged();
    m_bloomUp.reloadIfChanged();
    m_tonemap.reloadIfChanged();
}

std::string Renderer::shaderErrors() const {
    std::string all;
    for (const ShaderProgram* p : {&m_trace, &m_bloomDown, &m_bloomUp, &m_tonemap}) {
        if (!p->error().empty()) all += p->error();
    }
    return all;
}
