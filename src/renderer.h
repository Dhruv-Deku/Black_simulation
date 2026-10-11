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

// The central object as the ray tracer draws it (filled in by the app). Lengths and
// times in M.
struct CentralView {
    int kind = 0;                 // central::Kind
    double surface = 0.0;         // surface radius; 0 for a black hole
    double omega = 0.0;           // angular velocity of the surface
    double phase = 0.0;           // its rotation angle now (from phi_BL)
    double azimuthShift = 0.0;    // phi - phi_BL at the surface
    float temperature = 6000.0f;  // display temperature of the surface, K
    float hotSpotTemperature = 0.0f;  // magnetic polar caps (0: none)
    float brightness = 1.0f;
    glm::vec3 magneticAxis{0.0f, 0.0f, 1.0f};  // in the rotating frame
    bool pulsar = false, magnetar = false;
    bool detail = true;
    double coronaHeight = 0.0;    // scale height of a corona or atmosphere; 0 = none
    glm::vec3 coronaColor{0.0f};  // its color and strength (per unit length in M)
    double funnelRadius = 0.0;    // magnetospheric truncation radius; 0 = no accretion columns
    double diskInner = 6.0, diskOuter = 20.0;
    glm::vec3 sunDirection{0.6f, -0.5f, 0.62f};
    // Bodies that struck the surface: direction in the rotating frame and strength,
    // and how long ago it happened.
    std::vector<glm::vec4> impacts;
    std::vector<double> impactAges;
    double impactFade = 20.0;
    // A surface-wide flash (an X-ray burst, a nova, a collapse): age, strength, fade time.
    double flashAge = 1e30, flashStrength = 0.0, flashFade = 50.0;
};

struct FrameInputs {
    const Settings* settings = nullptr;
    const Camera* camera = nullptr;
    const GpuBodies* bodies = nullptr;
    const ParticleSystem* particles = nullptr;
    const CentralView* central = nullptr;
    int skipBody = -1;  // body the camera rides with, hidden
    double spin = 0.0;
    double scale = 1.0;  // render unit in M: the scene scale (see Settings)
    double time = 0.0;  // coordinate time, animates the disk
    int width = 1;      // output size in pixels
    int height = 1;
    float renderScale = 0.0f;  // overrides settings->renderScale when > 0
    float exposureBoost = 1.0f;  // a flash that lights up the whole view (a collapse)
};

// Ray traces the scene into an HDR buffer, adds bloom, tone maps it into an
// 8-bit image and copies that to the window. The GPU works in render units of
// scale M (central mass 1 / scale): everything is converted on the way in.
class Renderer {
public:
    bool init(const std::filesystem::path& shaderDir, std::string& error);
    void shutdown();

    // Sends new body-history samples (and any re-filled rows) to the GPU.
    void uploadHistory(History& history);
    // Replaces the lines drawn into the scene (points in M).
    void setPaths(const std::vector<PathLine>& lines, double spin, double scale);
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
    // The ray tracer for a black hole or a body at the center, with or without small
    // black holes among the bodies (shaders/trace.frag variants), built when first needed.
    ShaderProgram& traceProgram(bool centralBody, bool lenses);

    ShaderProgram m_trace[4];
    std::filesystem::path m_shaderDir;
    ShaderProgram m_bloomDown, m_bloomUp, m_tonemap;
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
    int m_historyGeneration = 0;
    GLuint m_pathSegments = 0, m_pathCells = 0, m_pathIndex = 0;  // SSBOs for the path lines
    GLuint m_pathOccupancy = 0;  // 3D texture: which path cells hold any segment
    int m_pathSegmentCount = 0;
    int m_skyWidth = 0;

    GLuint m_queries[kQueries] = {};
    int m_queryFrame = 0;
    float m_gpuMs = 0.0f;
};
