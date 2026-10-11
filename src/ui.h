#pragma once

#include "bodies.h"
#include "camera.h"
#include "particles.h"
#include "settings.h"
#include "trajectory.h"

#include <string>
#include <vector>

class Renderer;

// Camera paths for video recording.
enum RecordPath {
    RecordAsFlown = 0,  // whatever you do with the mouse and keys
    RecordOrbit = 1,    // one full turn around the hole
    RecordPushIn = 2,   // a slow approach toward the hole
};

struct VideoSize {
    const char* label;
    int width, height;
};
constexpr VideoSize kVideoSizes[] = {
    {"1280 x 720", 1280, 720}, {"1920 x 1080", 1920, 1080}, {"2560 x 1440", 2560, 1440}, {"3840 x 2160", 3840, 2160}};
constexpr int kVideoFps[] = {30, 60};

// Famous black holes, applied by App::applyScene.
enum Scene { SceneDefault, SceneM87, SceneSgrA, SceneGargantua, SceneCygnusX1, SceneCount };
constexpr const char* kSceneNames[] = {"Default", "M87*", "Sagittarius A*", "Gargantua", "Cygnus X-1"};
constexpr const char* kSceneDescriptions[] = {
    "The starting view.",
    "The first black hole ever imaged (Event Horizon Telescope, 2019), in the galaxy M87: 6.5 billion "
    "Suns, seen about 17 degrees from its spin axis, with a jet pointing almost at us.",
    "The black hole at the center of our Milky Way: 4 million Suns, a dim and turbulent disk.",
    "The black hole from the film Interstellar: nearly maximal spin, seen almost edge-on, with the "
    "Doppler effect switched off as in the film.",
    "The first black hole ever identified (1971): about 21 Suns, spinning near the maximum, "
    "feeding on a blue supergiant companion.",
};

// UI-only state that persists between frames.
struct UiState {
    bool showPanel = true;
    bool showOrbitMap = true;
    bool showHelp = false;

    // Custom launch form
    int launchKind = 0;  // BodyKind
    float launchR = 12.0f;
    float launchPolarDeg = 90.0f;
    float launchAzimuthDeg = 0.0f;
    float launchSpeed = 1.0f;  // multiple of the circular-orbit speed
    float launchTiltDeg = 0.0f;
    float launchRadial = 0.0f;
    bool launchPrograde = true;
    float launchRp = 8.0f;
    float launchRa = 20.0f;
    float launchTemperature = 5800.0f;
    float launchRadius = 0.5f;
    float launchMass = 1e-6f;  // units of the hole's mass

    // Orbit map
    float mapExtent = 30.0f;  // half-width of the view, in M
    int mapView = 0;          // 0 = from above, 1 = from the side

    // Video
    int videoSize = 1;
    int videoFps = 1;
    float videoQuality = 1.0f;
    float videoSeconds = 10.0f;
    int videoPath = RecordOrbit;

    char skyPath[512] = "";
    std::string skyStatus;
    std::string screenshotStatus;
    std::string videoStatus;
    std::string cameraStatus;

    // Launch planner: the predicted fate and path of the custom launch.
    bool plannerOpen = false;      // the launch panel is open this frame
    bool plannerValid = false;
    std::string plannerKey;        // the inputs the prediction was made for
    double plannerSpin = 0.0;
    kerr::Constants plannerConstants;
    kerr::State plannerState;
    trajectory::Analysis plannerAnalysis;
    trajectory::Path plannerPath;
    // With gravity between bodies: the planned body's path among the others (filled in
    // by App on a worker thread).
    bool plannerJointValid = false;
    trajectory::Path plannerJointPath;
    bool mapAiming = false;        // dragging out a launch on the orbit map

    // Ray inspector: what the light reaching a clicked pixel did.
    bool inspectorOpen = false;
    bool inspectValid = false;
    float inspectX = 0.0f, inspectY = 0.0f;  // clicked window position
    trajectory::RayReport inspect;
    std::string inspectBody;
};

struct UiContext {
    Settings& settings;
    Camera& camera;
    BodySystem& bodies;
    Renderer& renderer;
    ParticleSystem& particles;
    UiState& state;
    float fps = 0.0f;
    const char* gpuName = "";
    bool recording = false;
    int recordFrame = 0;
    int recordTotal = 0;

    // Predicted future paths of the bodies, indexed by Body::slot.
    const std::vector<trajectory::Path>* predictions = nullptr;
    double slowdown = 1.0;  // time runs this much slower (approaching the hole)

    // Requests carried out by the app after the UI is drawn.
    bool screenshotRequested = false;
    bool recordRequested = false;
    bool stopRecordingRequested = false;
    int cameraMotionRequested = -1;  // a CameraMotion
    int rideRequested = -1;          // history slot of a body to ride with
    int sceneRequested = -1;         // a Scene
    bool einsteinRingRequested = false;
};

// Turns the custom launch form into initial conditions. False if they are invalid.
bool plannedLaunch(const UiState& ui, double spin, kerr::Constants& c, kerr::State& s);

void drawUi(UiContext& ctx);
