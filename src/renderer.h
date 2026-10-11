#pragma once

#include "bodies.h"
#include "camera.h"
#include "gl_util.h"
#include "particles.h"
#include "settings.h"

#include <filesystem>
#include <string>
#include <vector>

// A polyline drawn into the scene as a thin glowing line (orbit trails, predictions).
struct PathLine {
    std::vector<glm::vec3> points;  // Cartesian
    glm::vec3 color{1.0f};
    float startOpacity = 1.0f;
    float endOpacity = 1.0f;
    bool dashed = false;
};

struct FrameInputs {
    const Settings* settings = nullptr;
    const Camera* camera = nullptr;
    const GpuBodies* bodies = nullptr;
    const ParticleSystem* particles = nullptr;
    int skipBody = -1;  // body the camera rides with, hidden
    double spin = 0.0;
    double time = 0.0;  // coordinate time, animates the disk
    int width = 1;      // output size in pixels
    int height = 1;
    float renderScale = 0.0f;  // overrides settings->renderScale when > 0
};

// Ray traces the scene into an HDR buffer, adds bloom, tone maps it into an
// 8-bit image and copies that to the window.
class Renderer {
public:
    bool init(const std::filesystem::path& shaderDir, std::string& error);
    void shutdown();

    // Sends new body-history samples (and any re-filled rows) to the GPU.
    void uploadHistory(History& history);
    // Replaces the lines drawn into the scene.
    void setPaths(const std::vector<PathLine>& lines, double spin);
    void render(const FrameInputs& in);
    // Shows the last image in the window, letterboxed if the aspect ratios differ.
    void present(int width, int height) const;
    void bindOutput() const { glBindFramebuffer(GL_FRAMEBUFFER, m_ldr.fbo); }
    bool saveScreenshot(const std::string& path, std::string& error) const;
    // The last image as top-down BGRA bytes (the layout video encoders take).
    void readOutput(std::vector<unsigned char>& bgra, int& width, int& height) const;

    bool loadSkyImage(const std::string& path, std::string& error);
    bool hasSkyImage() const { return m_skyTexture != 0; }

    void hotReload();
    std::string shaderErrors() const;
    float gpuMilliseconds() const { return m_gpuMs; }
    int traceWidth() const { return m_hdr.width; }
    int traceHeight() const { return m_hdr.height; }

private:
    static constexpr int kBloomLevels = 6;
    static constexpr int kQueries = 4;

    void drawFullscreen() const;
    void readTimer();

    ShaderProgram m_trace, m_bloomDown, m_bloomUp, m_tonemap;
    RenderTarget m_hdr, m_ldr;
    RenderTarget m_bloom[kBloomLevels];
    int m_bloomCount = 0;

    GLuint m_vao = 0;
    GLuint m_bodyUbo = 0;
    GLuint m_blackbodyTexture = 0;
    GLuint m_skyTexture = 0;
    GLuint m_placeholderTexture = 0;
    GLuint m_placeholderVolume = 0;
    GLuint m_placeholderOccupancy = 0;
    GLuint m_historyTexture = 0;
    long long m_historyUploaded = -1;
    GLuint m_pathSegments = 0, m_pathCells = 0, m_pathIndex = 0;  // SSBOs for the path lines
    GLuint m_pathOccupancy = 0;  // 3D texture: which path cells hold any segment
    int m_pathSegmentCount = 0;
    int m_skyWidth = 0;

    GLuint m_queries[kQueries] = {};
    int m_queryFrame = 0;
    float m_gpuMs = 0.0f;
};
