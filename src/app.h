#pragma once

#include "bodies.h"
#include "camera.h"
#include "particles.h"
#include "recorder.h"
#include "renderer.h"
#include "settings.h"
#include "ui.h"

#include <filesystem>
#include <future>
#include <optional>
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
    std::optional<int> relativity;
    std::optional<float> diskTemperature, diskBrightness;
    std::optional<float> debrisGlow, debrisThickness;
    std::optional<int> cameraMotion;  // CameraMotion
    std::optional<int> ride;          // index of a body to ride with
    std::optional<int> scene;         // Scene
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
    int debugView = 0;  // 1 = step heat map, 2 = escape direction

    // Video: record this many seconds to recordPath, then exit.
    std::string recordPath;
    double recordSeconds = 5.0;
    int recordFps = 30;
    int recordPathMode = 0;  // RecordPath

    bool showHelp = false;
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
    void setCameraMotion(CameraMotion motion);
    void rideWith(int slot);
    void syncRide(double dt);
    int ridingIndex() const;  // index in m_bodies.bodies of the ridden body, or -1
    // Time runs this much slower while the camera falls in close to the hole.
    double approachSlowdown() const;
    void updatePredictions();
    void buildPaths();
    void inspectPixel(double x, double y);
    void applyScene(int scene);
    void addEinsteinRing();

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
    } m_recording;

    std::vector<trajectory::Path> m_predictions = std::vector<trajectory::Path>(History::kSlots);  // by slot
    int m_predictNext = 0;
    // With gravity between bodies, all bodies (and the planned launch) are predicted
    // together on a worker thread; each result lands a moment later.
    std::future<std::vector<trajectory::Path>> m_jointJob;
    std::vector<int> m_jointSlots;   // slot of each path in the job (-1: the planned launch)
    std::string m_jointPlannerKey;   // planner inputs the job was started with

    float m_lastSpin = 0.0f;
    bool m_vsync = true;
    double m_scroll = 0.0;
    double m_mouseX = 0.0, m_mouseY = 0.0;
    bool m_dragLeft = false, m_dragRight = false;
    bool m_leftWas = false, m_middleWas = false;
    bool m_keyDown[512] = {};
    float m_fps = 0.0f;
    std::string m_gpuName;
};
