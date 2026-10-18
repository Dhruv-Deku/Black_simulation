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
    double turns = 0.0;           // turns it has made so far (counted on, modulo 4096)
    double azimuthShift = 0.0;    // phi - phi_BL at the surface
    float temperature = 6000.0f;  // display temperature of the surface, K
    float hotSpotTemperature = 0.0f;  // magnetic polar caps (0: none)
    float brightness = 1.0f;
    glm::vec3 magneticAxis{0.0f, 0.0f, 1.0f};  // in the rotating frame
    bool pulsar = false, magnetar = false;
    bool detail = true;
    // How the surface is drawn (surface::Look) and its scales (src/surface.h).
    int look = 0;
    double radiusKm = 0.0;        // the surface radius in km: the planet's air is drawn at true scale
    float planetSeed = 0.0f;      // which continents, storms, bands
    float planetLight = 1.0f;     // sunlight on a planet, relative to a star's disk
    bool atmosphere = false;      // an Earth-like sky
    float coronaScale = 0.0f;     // corona brightness per unit electron column (0: none)
    float chromosphereScale = 0.0f, chromosphereHeight = 0.002f;
    bool filmColor = false;       // stars in film colors (shaders/central_surface.glsl, filmHue)
    bool planetMaps = true;       // the Earth and Jupiter drawn from spacecraft maps, when loaded
    float envelope = 0.0f;        // a red giant's extended atmosphere, in film colors (0: none)
    float flowGlow = 0.0f;        // the thick turbulent inflow onto a neutron star, in film colors (0: none)
    bool sdoMap = false;          // a Sun-like star in film colors: its baked 30.4 nm map is needed
    float granuleFrequency = 0.0f, granuleContrast = 0.0f, granuleLane = 0.15f;
    float granulePhase = 0.0f;    // the convection cells' evolution now, in lifetimes
    double granuleRate = 0.0;     // lifetimes per M of coordinate time
    std::vector<glm::vec4> spots;      // direction in the rotating frame, angular radius
    std::vector<float> spotUmbra;      // umbra radius / spot radius
    float spotDrop = 0.0f;             // umbral temperature drop / surface temperature
    std::vector<glm::vec4> regions;    // active regions: direction, angular radius
    float capCos = 0.96f;              // cos of the polar caps' angular radius
    // White-light flares: direction in the rotating frame and strength, their ages, and
    // their light curve (rise and decay times, M).
    std::vector<glm::vec4> flares;
    std::vector<double> flareAges;
    double flareRise = 1.0, flareDecay = 10.0;
    double funnelRadius = 0.0;    // magnetospheric truncation radius; 0 = no accretion columns
    // A neutron star's magnetosphere (shaders/central_surface.glsl): its closed field
    // lines out to fieldReach (M), a magnetar's twist (radians between a line's
    // footpoints) and pair corona, and pulsar beams of half-opening angle beamWidth.
    float fieldGlow = 0.0f, fieldTwist = 0.0f, pairGlow = 0.0f;
    double fieldReach = 0.0;
    float beamGlow = 0.0f, beamWidth = 0.3f;
    // A cataclysmic variable: its companion (center in M, Roche-lobe radius in M; 0: none),
    // the direction to it, its temperature, the stream from L1 to the disk's rim (M) and
    // its half-width at L1 (M), the companion's azimuth and where the stream strikes.
    glm::vec4 companion{0.0f};
    glm::vec3 companionAxis{1.0f, 0.0f, 0.0f};
    float companionTemperature = 3100.0f;
    std::vector<glm::vec3> stream;
    double streamWidth = 0.0;
    float binaryPhase = 0.0f, spotPhase = 2.3f;
    double diskInner = 6.0, diskOuter = 20.0;
    glm::vec3 sunDirection{0.6f, -0.5f, 0.62f};
    // What lights the planets and moons among the bodies: 0, the disk (around a black
    // hole); 1, the star at the center, whose light at distance d (M) is
    // bodyLightScale / d^2 of the Sun's at the Earth; 2, the distant Sun along
    // sunDirection (around a planet), at bodyLightScale, in the planet's shadow blocked.
    int bodyLight = 0;
    float bodyLightTemperature = 5772.0f;  // display temperature of that light
    double bodyLightScale = 0.0;
    // Bodies that struck the surface: direction in the rotating frame and strength,
    // and how long ago it happened.
    std::vector<glm::vec4> impacts;
    std::vector<double> impactAges;
    double impactFade = 20.0;
    // A young neutron star drawn as films show one (central::Object::beacon): the glow
    // around it, and the starquakes ringing in its crust. For each quake: how long ago it
    // began (M), how far the magnetic poles swing at first (radians), the crust mode's
    // angular frequency (rad / M) and decay time (M); and which way it shakes the poles
    // (the star's rotating frame) with the mode's degree l.
    bool beacon = false;
    float haloGlow = 0.0f;
    std::vector<glm::vec4> quakes;
    std::vector<glm::vec4> quakeAxes;
    // A surface-wide flash (an X-ray burst, a nova, a collapse): age, strength, fade time.
    double flashAge = 1e30, flashStrength = 0.0, flashFade = 50.0;
};

// A tidal disruption's aftermath as the ray tracer draws it (filled in by App::flareView;
// see src/tde.h). Lengths and times in M.
struct FlareView {
    bool active = false;            // there is an outflow (or the light of one) to draw
    glm::vec3 axis{0.0f, 0.0f, 1.0f};  // its axis: the debris's orbital angular momentum (Cartesian)
    double launchRadius = 15.0;     // where the wind starts: the circularization radius
    double speed = 0.1;             // terminal speed on the axis (c)
    double extent = 0.0;            // how far out it is drawn (its front, capped)
    double clock = 0.0;             // coordinate time now, on the outflow's clock (from the disruption)
    double start = 0.0;             // on that clock, when the first debris came back
    double maxOpening = 1.0;        // the widest the outflow has been (radians)
    // How far from the axis gas (billows included) can be seen at radius R, for
    // kConeRows radii spaced evenly in ln R from launchRadius to extent: lets the ray
    // tracer skip at once the many steps that pass nowhere near it.
    static constexpr int kConeRows = 64;
    std::vector<float> cone;
    // History, row k at clock time tableStart + k tableStep: the rate at which mass is blown
    // out (relative to the fallback's peak), the opening half-angle (radians), the
    // escaping luminosity (relative to its peak), unused.
    double tableStart = 0.0, tableStep = 1.0;
    std::vector<glm::vec4> table;
    float thickness = 1.0f;         // opacity scale
    float brightness = 1.0f;        // how brightly the flare lights the outflow
    float falloff = 2.0f;           // the flare light falls off as R^-falloff
    float flareTemperature = 7000.0f;  // display temperature of the flare light (K)
    float glow = 1.0f;              // the hot base's own glow
    float envelope = 0.3f;          // the slow, quasi-spherical part of the wind
    float swirl = 0.0f;             // how far the wind turns about the axis on its way out (radians)
    // The small disk the returning gas forms: its brightness now (0: none), outer
    // radius and display temperature.
    float diskGlow = 0.0f;
    double diskOuter = 15.0;
    float diskTemperature = 6500.0f;
};

struct FrameInputs {
    const Settings* settings = nullptr;
    const Camera* camera = nullptr;
    const GpuBodies* bodies = nullptr;
    const ParticleSystem* particles = nullptr;
    const CentralView* central = nullptr;
    const FlareView* flare = nullptr;
    int skipBody = -1;  // body the camera rides with, hidden
    double spin = 0.0;
    double scale = 1.0;  // render unit in M: the scene scale (see Settings)
    double reach = 0.0;  // the farthest body from the center, in M: rays must reach it before escaping
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
    // farOrbits: around a star or planet the grid stretches out to the farthest orbit.
    void setPaths(const std::vector<PathLine>& lines, double spin, double scale, bool farOrbits = false);
    void render(const FrameInputs& in);
    // Shows the last image in the window, letterboxed if the aspect ratios differ.
    void present(int width, int height) const;
    void bindOutput() const { glBindFramebuffer(GL_FRAMEBUFFER, m_ldr.fbo); }
    bool saveScreenshot(const std::string& path, std::string& error) const;
    // The last image as top-down BGRA bytes (the layout video encoders take).
    void readOutput(std::vector<unsigned char>& bgra, int& width, int& height) const;

    bool loadSkyImage(const std::string& path, std::string& error);
    // Loads the Earth's and Jupiter's maps from assets/planets beside the shaders.
    void loadPlanetMaps();
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
    // The ray tracer for a black hole or a body at the center (one variant per kind of
    // surface), with or without small black holes or planets among the bodies
    // (shaders/trace.frag variants), built when first needed: code a scene does not use
    // would cost it registers, and so speed.
    ShaderProgram& traceProgram(int look, bool lenses, bool planets, bool outflow);

    static constexpr int kTraceLooks = 7;  // a black hole, and the six surface looks
    ShaderProgram m_trace[kTraceLooks * 8];
    GLuint m_outflowTexture = 0;  // FlareView::table (1D, RGBA32F)
    int m_outflowRows = 0;
    std::filesystem::path m_shaderDir;
    ShaderProgram m_bloomDown, m_bloomUp, m_tonemap;
    // The Sun's baked 30.4 nm map (shaders/sdo_bake.frag): a band of rows refreshed each
    // frame, all of it when it is first needed.
    ShaderProgram m_sdoBake;
    GLuint m_sdoMap = 0, m_sdoFbo = 0;
    int m_sdoBand = 0;
    bool m_sdoReady = false;
    static constexpr int kSdoWidth = 2048, kSdoHeight = 1024, kSdoBands = 16;
    void bakeSdo(const CentralView& cv, bool all);
    RenderTarget m_hdr, m_ldr;
    RenderTarget m_bloom[kBloomLevels];
    int m_bloomCount = 0;

    GLuint m_vao = 0;
    GLuint m_bodyUbo = 0;
    GLuint m_blackbodyTexture = 0;
    GLuint m_skyTexture = 0;
    // Spacecraft maps of the Earth (day, clouds, night lights) and Jupiter (assets/planets).
    GLuint m_earthDay = 0, m_earthClouds = 0, m_earthNight = 0, m_jupiterMap = 0;
    GLuint m_placeholderTexture = 0;
    GLuint m_placeholderVolume = 0;
    GLuint m_placeholderOccupancy = 0;
    GLuint m_historyTexture = 0;
    long long m_historyUploaded = -1;
    int m_historyGeneration = 0;
    GLuint m_pathSegments = 0, m_pathCells = 0, m_pathIndex = 0;  // SSBOs for the path lines
    GLuint m_pathOccupancy = 0;  // 3D texture: which path cells hold any segment
    double m_pathRMax = 100.0;  // outer radius of the path grid, render units
    int m_pathSegmentCount = 0;
    int m_skyWidth = 0;

    GLuint m_queries[kQueries] = {};
    int m_queryFrame = 0;
    float m_gpuMs = 0.0f;
};
