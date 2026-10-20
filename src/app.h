#pragma once

#include "bodies.h"
#include "camera.h"
#include "eht_view.h"
#include "json.h"
#include "merger_scene.h"
#include "particles.h"
#include "recorder.h"
#include "renderer.h"
#include "scene_sound.h"
#include "settings.h"
#include "timeline.h"
#include "ui.h"

#include <cstdint>
#include <filesystem>
#include <future>
#include <optional>
#include <random>
#include <string>
#include <vector>

struct GLFWwindow;

struct LaunchOptions {
    int width = 1600;
    int height = 900;
    std::string screenshotPath;  // if set: render, save this image and exit
    int frames = 4;              // frames rendered before that screenshot
    bool screenshotUi = false;   // include the interface in that screenshot
    double advance = 0.0;        // simulated time (M) to run before the first frame
    std::optional<float> spin;
    std::optional<double> distance, inclinationDeg, azimuthDeg;
    std::optional<float> fov;
    std::optional<int> sky;
    std::optional<float> renderScale;
    std::optional<int> upscaler;          // Settings::upscaler
    bool dynamicResolution = false;
    std::optional<int> relativity;
    std::optional<float> diskTemperature, diskBrightness;
    std::optional<float> debrisGlow, debrisThickness;
    std::optional<int> cameraMotion;  // CameraMotion
    std::optional<int> ride;          // index of a body to ride with
    std::optional<int> scene;         // Scene
    std::optional<int> central;       // central::catalog() entry
    std::optional<double> massSolar;  // physical mass of the central object
    std::optional<double> radiusKm;   // its radius
    std::optional<double> compactness;  // its radius in M (R / M)
    std::optional<double> surfaceTemperature, magneticField;
    bool physicalUnits = false;
    bool naked = false;  // a black hole may spin past 1: a naked singularity (--naked)
    std::optional<std::pair<double, double>> inspect;  // window position to inspect
    bool jets = false;
    bool noPaths = false;
    bool lookAhead = false;
    std::vector<std::string> presets;
    std::string skyImage;
    bool noDisk = false;
    bool noLightDelay = false;
    bool noGasPhysics = false;
    std::optional<float> gasCollisions, gasCooling;
    bool noBodyGravity = false;
    bool noFlare = false;
    bool starquake = false;  // a young neutron star's crust slips at the start
    int debugView = 0;  // 1 = step heat map, 2 = escape direction

    // The neutron star merger (src/merger_scene.h).
    bool merger = false;
    std::optional<int> mergerView;          // MergerScene::View
    std::optional<int> mergerParticles;     // per star
    std::optional<double> mergerCollapseMs;
    std::optional<double> mergerMass1, mergerMass2;
    double mergerAdvance = 0.0;             // seconds of viewing to run before the first frame
    bool noSound = false;
    bool noCaptions = false;
    bool sound = false;  // the scene's sound (--sound), live and in videos

    // Scenes (src/scene_file.h).
    std::string sceneFile;  // open this saved scene first (--scene-file)
    std::string saveScene;  // save the scene here once set up and advanced (--save-scene)
    // The timeline (src/timeline.h): on in the window, off for screenshots and videos unless
    // --seek goes back to a time (in the units of --advance) after the advance.
    std::optional<double> seek;
    bool noTimeline = false;
    // Camera keyframes (src/keyframes.h): hovering, looking at the center, at distance,
    // inclination and azimuth (degrees) as --distance and the rest take them, reached
    // after `seconds` from the one before.
    struct Keyframe {
        double distance, inclinationDeg, azimuthDeg, seconds;
    };
    std::vector<Keyframe> keyframes;
    bool keyframesPhysical = false;  // the camera moves along them for real (aberration)
    // What the Event Horizon Telescope sees (src/eht_view.h): EhtView::Target, or -1.
    int eht = -1;
    std::optional<double> ehtInclinationDeg;
    std::optional<int> ehtSize;
    std::optional<double> ehtHeight, ehtTemperature, ehtSpin;

    // Video: record this many seconds to recordPath, then exit.
    std::string recordPath;
    double recordSeconds = 5.0;
    int recordFps = 30;
    int recordPathMode = 0;  // RecordPath

    bool showHelp = false;
};

// What the window's callbacks leave for the app: scrolling, and files dropped onto it
// (paths in UTF-8).
struct WindowInbox {
    double scroll = 0.0;
    std::vector<std::string> dropped;
};

// Returns false with a message on bad arguments; sets showHelp for --help.
bool parseArguments(int argc, char** argv, LaunchOptions& out, std::string& error);
void printUsage();

class App {
public:
    int run(const LaunchOptions& options, const std::filesystem::path& shaderDir);

private:
    bool init(const LaunchOptions& options, const std::filesystem::path& shaderDir);
    void applyOptions(const LaunchOptions& options);
    void shutdown();
    void handleInput(float dt);
    bool pressed(int key);
    void saveScreenshot(const std::string& path);
    void advanceSimulation(double dt);
    // One frame of the simulation: everything that happens in simDt of coordinate time.
    // realDt (seconds) only eases the camera that follows a tidal disruption.
    void simulateFrame(double simDt, float realDt);

    // The timeline (src/timeline.h). Each frame it notices what the user changed, then
    // either steps the scene live by liveDt and records the step, or replays the past.
    void timelineFrame(float realDt, double liveDt);
    // A checkpoint of the scene as it is now (pinned: right after something the user did).
    void addCheckpoint(bool pinned);
    bool restoreCheckpoint(const Timeline::Checkpoint& checkpoint);
    // Restores and replays toward position `target` for at most `seconds` of wall time;
    // true once there.
    bool seekTowards(size_t target, double seconds);
    void releaseGasCopies(const std::vector<int>& slots);
    void clearTimeline();
    // Everything the user can change about the scene between frames, hashed: bodies, gas,
    // flashes, a moving camera. A change means they did something.
    std::uint64_t simulationFingerprint() const;
    UiContext::TimelineView timelineView() const;
    void timelineCommand(int command);
    // Camera keyframes (src/keyframes.h).
    void addKeyframe();
    // Moves the camera along the keyframes' path while it plays (or a video is filmed
    // along it); timeRate is the coordinate time per real second, for the camera's real
    // speed when it moves for real.
    void playKeyframes(float dt, double timeRate);
    void setCameraMotion(CameraMotion motion);
    void rideWith(int slot);
    void syncRide(double dt);
    int ridingIndex() const;  // index in m_bodies.bodies of the ridden body, or -1
    // Time runs this much slower while the camera falls in close to the hole.
    double approachSlowdown() const;
    // Time runs this much faster while a torn-apart star's debris is far out on its first
    // orbit, when little changes for hundreds of M.
    double flarePace() const;
    // Where the camera's distance goes while it follows a tidal disruption: close for the
    // star's approach and the disruption, back out with the stream and then the outflow.
    double flareCameraDistance() const;
    void updatePredictions();
    void buildPaths();
    void inspectPixel(double x, double y);
    void applyScene(int scene);
    void addEinsteinRing();
    // Puts a catalogue object at the center (camera, disk and scale follow).
    void applyCentral(int index);
    // Keeps everything in step with the central object: collapse, spin limits, scene
    // scale, what the bodies, gas and camera need to know, bursts.
    void syncCentral(double simDt, float realDt);
    CentralView centralView() const;
    // A tidal disruption's flare and outflow as the ray tracer draws them (src/tde.h).
    FlareView flareView() const;
    // A surface-wide flash (X-ray burst, nova, flare), seen with light delay.
    void triggerFlash(double strength, double fadeScene, const std::string& message);
    double sceneScale() const { return m_settings.central.sceneScale(); }

    // Scenes (src/scene_file.h): everything that decides what happens next, in two files.
    // A loaded scene carries on exactly as the saved one would have.
    bool saveScene(const std::filesystem::path& path);
    bool loadScene(const std::filesystem::path& path);
    // Ctrl+S saves to the scene last saved or opened (choose: always ask where), Ctrl+O opens.
    void saveSceneInteractive(bool choose);
    void openSceneInteractive();
    // The app's own share of a scene: flashes, flares, starquakes, its random numbers.
    json::Value saveAppState() const;
    bool loadAppState(const json::Value& state, bool checkOnly, std::string& error);

    bool startRecording(const std::string& path, int width, int height, int fps, float quality, int cameraPath,
                        double seconds);
    void recordFrame();
    void stopRecording();

    GLFWwindow* m_window = nullptr;
    Settings m_settings;
    Camera m_camera;
    BodySystem m_bodies;
    ParticleSystem m_particles;
    Renderer m_renderer;
    UiState m_ui;
    MergerScene m_merger;

    struct Recording {
        VideoRecorder encoder;
        bool active = false;
        std::string path;
        int width = 0, height = 0, fps = 30;
        float quality = 1.0f;
        int cameraPath = 0;
        int frame = 0, total = 0;
        double startPhi = 0.0, startR = 0.0;
        std::vector<unsigned char> pixels;
        double soundSquares = 0.0;  // sum of the sound track's squared samples, and how many
        long long soundSamples = 0;
    } m_recording;

    std::vector<trajectory::Path> m_predictions = std::vector<trajectory::Path>(History::kSlots);  // by slot
    int m_predictNext = 0;
    // With gravity between bodies, all bodies (and the planned launch) are predicted
    // together on a worker thread; each result lands a moment later.
    std::future<std::vector<trajectory::Path>> m_jointJob;
    std::vector<int> m_jointSlots;   // slot of each path in the job (-1: the planned launch)
    std::string m_jointPlannerKey;   // planner inputs the job was started with

    float m_lastSpin = 0.0f;
    double m_lastScale = 1.0;
    double m_flashTime = -1e30, m_flashStrength = 0.0, m_flashFade = 50.0;  // coordinate times (M)
    double m_burstFuel = 0.0;        // fraction of the way to the next thermonuclear flash
    // White-light flares on an active star: where (rotating frame), when (M), how bright.
    struct StarFlare {
        glm::vec3 direction;
        double time;
        float strength;
    };
    std::vector<StarFlare> m_flares;
    // Starquakes on a young neutron star (central::Object::beacon): when, how far the
    // magnetic poles swing at first (radians), the crust mode's degree, and which way it
    // shakes them (the star's rotating frame, across the magnetic axis).
    struct Quake {
        double time;
        float amplitude;
        int l;
        glm::vec3 direction;
    };
    std::vector<Quake> m_quakes;
    // The crust slips now: a glitch, and the crust rings.
    void starquake();
    double m_lastAccreted = 0.0;     // bodies and gas that had hit the surface at the last check
    float m_screenFlash = 0.0f;      // a collapse lights up the whole view for a moment
    // The stream of a cataclysmic variable, cached (centralView): mass ratio, rim, path.
    mutable double m_streamQ = -1.0, m_streamRim = -1.0;
    mutable std::vector<std::array<double, 2>> m_streamPath;
    std::mt19937 m_rng{4242};
    bool m_vsync = true;
    WindowInbox m_inbox;
    keyframes::Track m_track;
    EhtView m_eht;
    SceneSound m_sound;
    float m_dynamicScale = 0.0f;  // the resolution dynamic resolution is heading for (smoothed)
    // What each body, pulsar and heavy body sounds like now (src/scene_sound.h), and the
    // rumble of the tide at the camera. timeRate: coordinate time per real second.
    std::vector<SceneSound::Voice> soundVoices(double timeRate, double& rumble, double& waveHz) const;
    double m_waveHz = 0.0;  // the true frequency of the loudest gravitational wave (Hz), for the panel
    int m_ehtSizeOverride = 0;  // --eht-size
    // The telescope view's options from the panel, the central object and the camera.
    EhtView::Options ehtOptions() const;
    Timeline m_timeline;
    bool m_timelineOn = true;
    std::vector<ParticleSystem::GpuCopy> m_gasCopies;  // by Timeline::Checkpoint::gpu (buffer 0: a free slot)
    std::uint64_t m_fingerprint = 0;  // of the scene at the end of the last frame
    long long m_seekTarget = -1;      // a position being replayed toward (-1: none)
    bool m_reverse = false;           // playing the past backwards
    std::filesystem::path m_scenePath;  // the scene last saved or opened
    double m_mouseX = 0.0, m_mouseY = 0.0;
    bool m_dragLeft = false, m_dragRight = false;
    bool m_leftWas = false, m_middleWas = false;
    bool m_keyDown[512] = {};
    float m_fps = 0.0f;
    std::string m_gpuName;
};
