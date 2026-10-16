#include "app.h"

#include "gravity.h"
#include "kerr.h"

#include <glad/gl.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <map>

namespace {

const std::map<std::string, Preset>& presetNames() {
    static const std::map<std::string, Preset> names = {
        {"isco", Preset::IscoOrbit},        {"circular", Preset::StableOrbit},  {"retrograde", Preset::RetrogradeOrbit},
        {"precessing", Preset::Precessing}, {"zoom-whirl", Preset::ZoomWhirl},  {"polar", Preset::PolarOrbit},
        {"plunge", Preset::Plunge},         {"flyby", Preset::Flyby},           {"cluster", Preset::Cluster},
        {"planet", Preset::RockyPlanet},    {"gas-giant", Preset::GasGiant},    {"pulsar", Preset::Pulsar},
        {"tde", Preset::TidalDisruption},   {"cloud", Preset::GasCloud},        {"ring", Preset::DebrisRing},
        {"grazing", Preset::GrazingStar},   {"binary", Preset::BinaryStar},     {"moon", Preset::PlanetWithMoon},
    };
    return names;
}

void APIENTRY glDebugCallback(GLenum, GLenum type, GLuint, GLenum severity, GLsizei, const GLchar* message,
                              const void*) {
    if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return;
    if (type == GL_DEBUG_TYPE_ERROR || severity == GL_DEBUG_SEVERITY_HIGH) {
        std::fprintf(stderr, "[GL] %s\n", message);
    }
}

void scrollCallback(GLFWwindow* window, double, double dy) {
    if (auto* scroll = static_cast<double*>(glfwGetWindowUserPointer(window))) *scroll += dy;
}

std::string timestampedPath(const char* folder, const char* extension) {
    std::time_t now = std::time(nullptr);
    std::tm local{};
#ifdef _WIN32
    localtime_s(&local, &now);
#else
    localtime_r(&now, &local);
#endif
    char stamp[32];
    std::strftime(stamp, sizeof(stamp), "%Y%m%d-%H%M%S", &local);
    return std::string(folder) + "/blackhole-" + stamp + extension;
}

void ensureParent(const std::string& path) {
    std::error_code ec;
    std::filesystem::path parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent, ec);
}

void applyStyle(float scale) {
    ImGuiStyle& style = ImGui::GetStyle();
    ImGui::StyleColorsDark();
    style.WindowRounding = 6.0f;
    style.FrameRounding = 4.0f;
    style.GrabRounding = 4.0f;
    style.WindowBorderSize = 0.0f;
    ImVec4* c = style.Colors;
    c[ImGuiCol_WindowBg] = ImVec4(0.06f, 0.065f, 0.08f, 0.88f);
    c[ImGuiCol_Header] = ImVec4(0.55f, 0.32f, 0.12f, 0.45f);
    c[ImGuiCol_HeaderHovered] = ImVec4(0.75f, 0.45f, 0.15f, 0.6f);
    c[ImGuiCol_HeaderActive] = ImVec4(0.85f, 0.5f, 0.18f, 0.7f);
    c[ImGuiCol_FrameBg] = ImVec4(0.14f, 0.15f, 0.19f, 0.9f);
    c[ImGuiCol_SliderGrab] = ImVec4(0.95f, 0.6f, 0.25f, 0.9f);
    c[ImGuiCol_SliderGrabActive] = ImVec4(1.0f, 0.72f, 0.35f, 1.0f);
    c[ImGuiCol_CheckMark] = ImVec4(1.0f, 0.68f, 0.3f, 1.0f);
    c[ImGuiCol_Button] = ImVec4(0.55f, 0.32f, 0.12f, 0.5f);
    c[ImGuiCol_ButtonHovered] = ImVec4(0.75f, 0.45f, 0.15f, 0.7f);
    c[ImGuiCol_ButtonActive] = ImVec4(0.85f, 0.5f, 0.18f, 0.9f);
    c[ImGuiCol_PlotHistogram] = ImVec4(0.95f, 0.6f, 0.25f, 0.9f);
    c[ImGuiCol_TitleBgActive] = ImVec4(0.12f, 0.1f, 0.09f, 1.0f);
    style.ScaleAllSizes(scale);
    style.FontScaleDpi = scale;
}

}  // namespace

// --- Command line -------------------------------------------------------------------

void printUsage() {
    std::printf(
        "Black Hole Simulator: real-time ray tracing of a rotating (Kerr) black hole.\n\n"
        "usage: blackhole [options]\n"
        "  --size WxH             window size (default 1600x900)\n"
        "  --spin A               black hole spin a/M in [0, 0.999]\n"
        "  --distance R           camera distance in M\n"
        "  --inclination DEG      camera angle from the spin axis (90 = edge-on)\n"
        "  --azimuth DEG          camera azimuth\n"
        "  --fov DEG              vertical field of view\n"
        "  --camera hover|orbit|fall\n"
        "                         hold the camera still, put it in orbit, or drop it\n"
        "  --ride N               ride along with body number N (0 = first added)\n"
        "  --look-ahead           while moving, face the direction of travel\n"
        "  --scene default|m87|sgra|gargantua|cygnus\n"
        "                         famous black holes\n"
        "  --jets                 relativistic jets along the spin axis\n"
        "  --no-paths             do not draw orbit trails and predictions in the view\n"
        "  --inspect X,Y          trace the light reaching window pixel X,Y and print what it did\n"
        "  --sky stars|grid|image|black\n"
        "  --sky-image PATH       equirectangular background image (.jpg, .png, .hdr)\n"
        "  --relativity full|beaming|off\n"
        "  --no-light-delay       draw things where they are now, not where their light left them\n"
        "  --no-gas-physics       debris without pressure, shocks, heating or cooling\n"
        "  --gas-collisions K     how readily crossing gas streams shock (default 0.3)\n"
        "  --gas-cooling T        gas cooling time in M at 10,000 K (default 60)\n"
        "  --no-gravity           bodies do not pull on each other or on the gas\n"
        "  --disk-temp K          peak disk temperature in kelvin\n"
        "  --disk-brightness B    disk brightness multiplier\n"
        "  --debris-glow B        debris brightness multiplier\n"
        "  --debris-thickness K   debris opacity (low = translucent glow)\n"
        "  --preset NAME          add objects (may repeat):\n"
        "                           orbits:  isco circular retrograde precessing zoom-whirl\n"
        "                                    polar plunge flyby cluster\n"
        "                           objects: planet gas-giant pulsar tde grazing binary moon\n"
        "                                    cloud ring einstein-ring\n"
        "  --advance T            run the simulation for T M before starting\n"
        "  --scale S              ray tracing resolution scale (0.25 to 2)\n"
        "  --no-disk              hide the accretion disk\n"
        "  --heatmap              show integration steps per pixel instead of the image\n"
        "  --debug-directions     show each ray's final direction as a color\n"
        "  --screenshot PATH      render without the UI, save a PNG and exit\n"
        "  --frames N             frames to render before that screenshot (default 4)\n"
        "  --ui                   include the interface in that screenshot\n"
        "  --record PATH          record an MP4 video (at --size) and exit\n"
        "  --duration S           video length in seconds (default 5)\n"
        "  --fps N                video frame rate (default 30)\n"
        "  --camera-path flown|orbit|push-in\n"
        "                         how the camera moves during the video\n");
}

bool parseArguments(int argc, char** argv, LaunchOptions& out, std::string& error) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto value = [&]() -> const char* {
            if (i + 1 >= argc) {
                error = arg + " needs a value";
                return nullptr;
            }
            return argv[++i];
        };
        auto number = [&](double& v) {
            const char* s = value();
            if (!s) return false;
            char* end = nullptr;
            v = std::strtod(s, &end);
            if (end == s || *end != '\0') {
                error = "bad number for " + arg + ": " + s;
                return false;
            }
            return true;
        };
        // Matches a word against a list; the index of the match, or -1.
        auto choice = [&](std::initializer_list<const char*> words) {
            const char* s = value();
            if (!s) return -1;
            int index = 0;
            for (const char* w : words) {
                if (std::string(s) == w) return index;
                ++index;
            }
            error = "unknown value for " + arg + ": " + s;
            return -1;
        };
        double v = 0.0;
        if (arg == "--help" || arg == "-h") {
            out.showHelp = true;
        } else if (arg == "--size") {
            const char* s = value();
            if (!s || std::sscanf(s, "%dx%d", &out.width, &out.height) != 2 || out.width < 16 || out.height < 16) {
                if (error.empty()) error = "--size expects WxH, e.g. 1920x1080";
                return false;
            }
        } else if (arg == "--spin") {
            if (!number(v)) return false;
            out.spin = static_cast<float>(std::clamp(v, 0.0, 0.999));
        } else if (arg == "--distance") {
            if (!number(v)) return false;
            out.distance = v;
        } else if (arg == "--inclination") {
            if (!number(v)) return false;
            out.inclinationDeg = v;
        } else if (arg == "--azimuth") {
            if (!number(v)) return false;
            out.azimuthDeg = v;
        } else if (arg == "--fov") {
            if (!number(v)) return false;
            out.fov = static_cast<float>(v);
        } else if (arg == "--scale") {
            if (!number(v)) return false;
            out.renderScale = static_cast<float>(std::clamp(v, 0.1, 2.0));
        } else if (arg == "--disk-temp") {
            if (!number(v)) return false;
            out.diskTemperature = static_cast<float>(std::clamp(v, 1000.0, 100000.0));
        } else if (arg == "--disk-brightness") {
            if (!number(v)) return false;
            out.diskBrightness = static_cast<float>(std::max(v, 0.0));
        } else if (arg == "--debris-glow") {
            if (!number(v)) return false;
            out.debrisGlow = static_cast<float>(std::max(v, 0.0));
        } else if (arg == "--debris-thickness") {
            if (!number(v)) return false;
            out.debrisThickness = static_cast<float>(std::max(v, 0.0));
        } else if (arg == "--advance") {
            if (!number(v)) return false;
            out.advance = std::max(v, 0.0);
        } else if (arg == "--frames") {
            if (!number(v)) return false;
            out.frames = std::max(1, static_cast<int>(v));
        } else if (arg == "--duration") {
            if (!number(v)) return false;
            out.recordSeconds = std::clamp(v, 0.1, 600.0);
        } else if (arg == "--fps") {
            if (!number(v)) return false;
            out.recordFps = std::clamp(static_cast<int>(v), 1, 120);
        } else if (arg == "--sky") {
            int c = choice({"stars", "grid", "image", "black"});
            if (c < 0) return false;
            out.sky = c;
        } else if (arg == "--relativity") {
            int c = choice({"full", "beaming", "off"});
            if (c < 0) return false;
            out.relativity = c;
        } else if (arg == "--camera") {
            int c = choice({"hover", "orbit", "fall"});
            if (c < 0) return false;
            out.cameraMotion = c;
        } else if (arg == "--camera-path") {
            int c = choice({"flown", "orbit", "push-in"});
            if (c < 0) return false;
            out.recordPathMode = c;
        } else if (arg == "--scene") {
            int c = choice({"default", "m87", "sgra", "gargantua", "cygnus"});
            if (c < 0) return false;
            out.scene = c;
        } else if (arg == "--ride") {
            if (!number(v)) return false;
            out.ride = static_cast<int>(v);
        } else if (arg == "--inspect") {
            const char* s = value();
            double x = 0.0, y = 0.0;
            if (!s || std::sscanf(s, "%lf,%lf", &x, &y) != 2) {
                if (error.empty()) error = "--inspect expects X,Y (window pixels)";
                return false;
            }
            out.inspect = std::make_pair(x, y);
        } else if (arg == "--jets") {
            out.jets = true;
        } else if (arg == "--no-paths") {
            out.noPaths = true;
        } else if (arg == "--look-ahead") {
            out.lookAhead = true;
        } else if (arg == "--sky-image") {
            const char* s = value();
            if (!s) return false;
            out.skyImage = s;
        } else if (arg == "--preset") {
            const char* s = value();
            if (!s) return false;
            if (!presetNames().count(s) && std::string(s) != "einstein-ring") {
                error = std::string("unknown preset: ") + s;
                return false;
            }
            out.presets.push_back(s);
        } else if (arg == "--screenshot") {
            const char* s = value();
            if (!s) return false;
            out.screenshotPath = s;
        } else if (arg == "--record") {
            const char* s = value();
            if (!s) return false;
            out.recordPath = s;
        } else if (arg == "--no-disk") {
            out.noDisk = true;
        } else if (arg == "--no-light-delay") {
            out.noLightDelay = true;
        } else if (arg == "--no-gas-physics") {
            out.noGasPhysics = true;
        } else if (arg == "--gas-collisions") {
            if (!number(v)) return false;
            out.gasCollisions = static_cast<float>(std::max(v, 0.0));
        } else if (arg == "--gas-cooling") {
            if (!number(v)) return false;
            out.gasCooling = static_cast<float>(std::max(v, 1.0));
        } else if (arg == "--no-gravity") {
            out.noBodyGravity = true;
        } else if (arg == "--heatmap") {
            out.debugView = 1;
        } else if (arg == "--debug-directions") {
            out.debugView = 2;
        } else if (arg == "--ui") {
            out.screenshotUi = true;
        } else {
            error = "unknown option: " + arg;
            return false;
        }
    }
    return true;
}

// --- App ---------------------------------------------------------------------------------

bool App::init(const LaunchOptions& options, const std::filesystem::path& shaderDir) {
    glfwSetErrorCallback([](int code, const char* text) { std::fprintf(stderr, "[GLFW %d] %s\n", code, text); });
    if (!glfwInit()) return false;

    const bool headless = !options.screenshotPath.empty() || !options.recordPath.empty();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Screenshots need the window to be exactly the requested size.
    glfwWindowHint(GLFW_SCALE_TO_MONITOR, headless ? GLFW_FALSE : GLFW_TRUE);
    glfwWindowHint(GLFW_VISIBLE, headless ? GLFW_FALSE : GLFW_TRUE);
#ifndef NDEBUG
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
#endif
    m_window = glfwCreateWindow(options.width, options.height, "Black Hole Simulator", nullptr, nullptr);
    if (!m_window) {
        std::fprintf(stderr, "Could not create an OpenGL 4.6 window. Is the graphics driver up to date?\n");
        return false;
    }
    glfwMakeContextCurrent(m_window);
    if (!gladLoadGL(glfwGetProcAddress)) {
        std::fprintf(stderr, "Could not load OpenGL functions.\n");
        return false;
    }
    glfwSwapInterval(headless ? 0 : 1);

    m_gpuName = reinterpret_cast<const char*>(glGetString(GL_RENDERER));
    std::printf("OpenGL %s on %s\n", glGetString(GL_VERSION), m_gpuName.c_str());
    glEnable(GL_DEBUG_OUTPUT);
    glDebugMessageCallback(glDebugCallback, nullptr);

    // Our scroll handler must be installed before ImGui's, which chains to it.
    glfwSetWindowUserPointer(m_window, &m_scroll);
    glfwSetScrollCallback(m_window, scrollCallback);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    if (headless) io.IniFilename = nullptr;
    float xscale = 1.0f, yscale = 1.0f;
    glfwGetWindowContentScale(m_window, &xscale, &yscale);
    applyStyle(std::max(xscale, 1.0f));
#ifdef _WIN32
    // ImGui's built-in pixel font looks rough when scaled; prefer Segoe UI.
    const char* font = "C:\\Windows\\Fonts\\segoeui.ttf";
    if (std::filesystem::exists(font)) io.Fonts->AddFontFromFileTTF(font, 16.0f);
#endif
    ImGui_ImplGlfw_InitForOpenGL(m_window, true);
    ImGui_ImplOpenGL3_Init("#version 460");

    std::string error;
    std::printf("Shaders: %s\n", shaderDir.string().c_str());
    if (!m_renderer.init(shaderDir, error)) {
        std::fprintf(stderr, "Shader errors (fix and save; they reload automatically):\n%s\n", error.c_str());
    }
    if (!m_particles.init(shaderDir, error)) {
        std::fprintf(stderr, "Debris shader errors:\n%s\n", error.c_str());
    }
    applyOptions(options);
    return true;
}

void App::applyOptions(const LaunchOptions& o) {
    if (o.scene) applyScene(*o.scene);
    if (o.spin) m_settings.spin = *o.spin;
    m_lastSpin = m_settings.spin;
    if (o.distance) m_camera.r = *o.distance;
    if (o.inclinationDeg) m_camera.theta = std::clamp(*o.inclinationDeg, 1.0, 179.0) * kerr::kPi / 180.0;
    if (o.azimuthDeg) m_camera.phi = *o.azimuthDeg * kerr::kPi / 180.0;
    if (o.fov) m_camera.fovDeg = *o.fov;
    m_camera.r = std::max(m_camera.r, kerr::horizonRadius(m_settings.spin) + 0.3);
    if (o.sky) m_settings.skyMode = *o.sky;
    if (o.renderScale) m_settings.renderScale = *o.renderScale;
    if (o.relativity) m_settings.relativity = *o.relativity;
    if (o.diskTemperature) m_settings.diskTemperature = *o.diskTemperature;
    if (o.diskBrightness) m_settings.diskBrightness = *o.diskBrightness;
    if (o.debrisGlow) m_settings.debrisBrightness = *o.debrisGlow;
    if (o.debrisThickness) m_settings.debrisOpacity = *o.debrisThickness;
    if (o.noDisk) m_settings.diskEnabled = false;
    if (o.noLightDelay) m_settings.lightDelay = false;
    if (o.noGasPhysics) m_settings.gasPhysics = false;
    if (o.gasCollisions) m_settings.gasCollisions = *o.gasCollisions;
    if (o.gasCooling) m_settings.gasCoolingTime = *o.gasCooling;
    if (o.noBodyGravity) m_settings.bodyGravity = false;
    if (o.jets) m_settings.jetsEnabled = true;
    if (o.noPaths) m_settings.showPaths = false;
    if (o.debugView) m_settings.debugView = o.debugView;
    if (!o.skyImage.empty()) {
        std::string error;
        if (m_renderer.loadSkyImage(o.skyImage, error)) {
            if (!o.sky) m_settings.skyMode = 2;
        } else {
            std::fprintf(stderr, "Could not load %s: %s\n", o.skyImage.c_str(), error.c_str());
        }
    }
    for (const std::string& name : o.presets) {
        if (name == "einstein-ring") addEinsteinRing();
        else m_bodies.addPreset(m_settings.spin, presetNames().at(name));
    }
    if (o.cameraMotion) setCameraMotion(static_cast<CameraMotion>(*o.cameraMotion));
    if (o.ride && *o.ride >= 0 && *o.ride < static_cast<int>(m_bodies.bodies.size())) {
        rideWith(m_bodies.bodies[*o.ride].slot);
    }
    for (double t = 0.0; t < o.advance; t += 1.0) {
        const double dt = std::min(1.0, o.advance - t);
        advanceSimulation(dt);
        syncRide(dt);
    }
    if (o.lookAhead) {
        m_camera.lookAhead = true;
        if (m_camera.motion != CameraMotion::Hover) m_camera.faceVelocity(m_settings.spin);
    }
    if (o.inspect) {
        inspectPixel(o.inspect->first, o.inspect->second);
        const trajectory::RayReport& r = m_ui.inspect;
        std::printf("Inspected pixel (%.0f, %.0f): light came from %s\n", o.inspect->first, o.inspect->second,
                    trajectory::rayEndName(r.end));
        std::printf("  r_end = %.3f  g = %.3f  turns = %.2f  r_min = %.3f  delay = %.1f M  b = %.3f  q = %.3f\n",
                    r.rEnd, r.g, r.turns, r.rMin, r.delay, r.b, r.q);
    }
}

void App::shutdown() {
    if (m_recording.active) stopRecording();
    m_particles.shutdown();
    m_renderer.shutdown();
    if (ImGui::GetCurrentContext()) {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
    }
    if (m_window) glfwDestroyWindow(m_window);
    glfwTerminate();
}

bool App::pressed(int key) {
    bool down = glfwGetKey(m_window, key) == GLFW_PRESS;
    bool was = m_keyDown[key];
    m_keyDown[key] = down;
    return down && !was;
}

void App::handleInput(float dt) {
    ImGuiIO& io = ImGui::GetIO();
    const double minR = kerr::horizonRadius(m_settings.spin) + 0.3;
    const bool moving = m_camera.motion != CameraMotion::Hover;

    // Mouse: a drag belongs to the camera only if it started outside the UI.
    double x, y;
    glfwGetCursorPos(m_window, &x, &y);
    const double dx = x - m_mouseX, dy = y - m_mouseY;
    m_mouseX = x;
    m_mouseY = y;
    const bool left = glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    const bool right = glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    const bool middle = glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;
    const bool ctrl = glfwGetKey(m_window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
                      glfwGetKey(m_window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
    // Ctrl+click or middle-click: inspect the light reaching that pixel (when the ray
    // inspector is on; otherwise Ctrl+drag orbits like a plain drag).
    const bool inspecting = m_settings.rayInspector && ctrl;
    if (m_settings.rayInspector && !io.WantCaptureMouse &&
        ((middle && !m_middleWas) || (left && !m_leftWas && ctrl))) {
        inspectPixel(x, y);
    }
    m_middleWas = middle;
    m_leftWas = left;
    if (!left) m_dragLeft = false;
    else if (!m_dragLeft && !io.WantCaptureMouse && !inspecting) m_dragLeft = true;
    if (!right) m_dragRight = false;
    else if (!m_dragRight && !io.WantCaptureMouse) m_dragRight = true;

    const double sensitivity = 0.005 * std::clamp(m_camera.fovDeg / 60.0, 0.3, 2.0);
    // While riding an orbit or falling, the path sets the position: dragging looks around.
    if (m_dragLeft && !moving) m_camera.orbit(-dx * sensitivity, -dy * sensitivity);
    if (m_dragRight || (m_dragLeft && moving)) {
        m_camera.look(static_cast<float>(dx * sensitivity), static_cast<float>(-dy * sensitivity));
    }
    if (!io.WantCaptureMouse && m_scroll != 0.0 && !moving) m_camera.zoom(std::pow(0.9, m_scroll), minR);
    m_scroll = 0.0;

    if (io.WantCaptureKeyboard) return;
    if (pressed(GLFW_KEY_SPACE)) m_settings.paused = !m_settings.paused;
    if (pressed(GLFW_KEY_TAB)) m_ui.showPanel = !m_ui.showPanel;
    if (pressed(GLFW_KEY_H)) m_ui.showHelp = !m_ui.showHelp;
    if (pressed(GLFW_KEY_R)) m_camera.reset();
    if (pressed(GLFW_KEY_G)) m_settings.skyMode = m_settings.skyMode == 1 ? 0 : 1;
    if (pressed(GLFW_KEY_I)) m_settings.rayInspector = !m_settings.rayInspector;
    if (pressed(GLFW_KEY_F12)) saveScreenshot(timestampedPath("screenshots", ".png"));
    if (pressed(GLFW_KEY_ESCAPE) && m_recording.active) stopRecording();

    if (moving) return;
    auto held = [&](int key) { return glfwGetKey(m_window, key) == GLFW_PRESS; };
    const double speed = 1.2 * dt;
    if (held(GLFW_KEY_W)) m_camera.zoom(std::exp(-speed), minR);
    if (held(GLFW_KEY_S)) m_camera.zoom(std::exp(speed), minR);
    if (held(GLFW_KEY_A)) m_camera.orbit(speed, 0.0);
    if (held(GLFW_KEY_D)) m_camera.orbit(-speed, 0.0);
    if (held(GLFW_KEY_Q)) m_camera.orbit(0.0, -speed);
    if (held(GLFW_KEY_E)) m_camera.orbit(0.0, speed);
}

void App::saveScreenshot(const std::string& path) {
    ensureParent(path);
    std::string error;
    if (m_renderer.saveScreenshot(path, error)) {
        m_ui.screenshotStatus = "Saved " + path;
        std::printf("Saved %s\n", path.c_str());
    } else {
        m_ui.screenshotStatus = "Failed: " + error;
        std::fprintf(stderr, "Screenshot failed: %s\n", error.c_str());
    }
}

void App::setCameraMotion(CameraMotion motion) {
    std::string message;
    switch (motion) {
        case CameraMotion::Hover:
            m_camera.hover(m_settings.spin, message);
            break;
        case CameraMotion::Orbit:
            m_camera.startOrbit(m_settings.spin, message);
            break;
        case CameraMotion::FreeFall:
            m_camera.startFall(m_settings.spin, message);
            break;
    }
    m_ui.cameraStatus = message;
}

void App::rideWith(int slot) {
    for (const Body& b : m_bodies.bodies) {
        if (b.slot != slot) continue;
        m_camera.rideWith(slot, b.constants, b.state, m_settings.spin, 0.0);
        m_camera.yaw = m_camera.pitch = 0.0f;
        m_ui.cameraStatus = "Riding with " + b.name + ".";
        return;
    }
}

double App::approachSlowdown() const {
    // Falling in, the last few M before the horizon and the whole way down to the inner
    // horizon take only a few M of time: slow down so the crossing can be watched.
    if (!m_settings.slowNearHole || m_camera.ended) return 1.0;
    const bool falling = m_camera.motion == CameraMotion::FreeFall ||
                         (m_camera.motion == CameraMotion::RideAlong && m_camera.path.vr < 0.0);
    if (!falling) return 1.0;
    const double x = m_camera.r / (4.0 * kerr::horizonRadius(m_settings.spin));
    return std::clamp(x * x, 0.04, 1.0);
}

int App::ridingIndex() const {
    if (m_camera.motion != CameraMotion::RideAlong) return -1;
    for (size_t i = 0; i < m_bodies.bodies.size(); ++i) {
        if (m_bodies.bodies[i].slot == m_camera.rideSlot) return static_cast<int>(i);
    }
    return -1;
}

void App::syncRide(double dt) {
    if (m_camera.motion != CameraMotion::RideAlong) return;
    const int i = ridingIndex();
    if (i < 0) {
        std::string message;
        if (m_camera.hover(m_settings.spin, message)) {
            m_ui.cameraStatus = "The body you were riding is gone; hovering where it was.";
        } else {
            m_camera.coast(m_camera.geodesic, m_camera.path);
            m_ui.cameraStatus = "The body you were riding is gone; you fall on along its path.";
        }
        return;
    }
    const Body& b = m_bodies.bodies[i];
    m_camera.rideWith(b.slot, b.constants, b.state, m_settings.spin, dt);
}

void App::updatePredictions() {
    const auto& bodies = m_bodies.bodies;
    const double spin = m_settings.spin;
    std::vector<gravity::Member> members;
    std::vector<int> slots;
    for (const Body& b : bodies) {
        members.push_back({b.constants, b.state, b.mass, 0.1 * b.radius});
        slots.push_back(b.slot);
    }

    if (m_settings.bodyGravity && gravity::interacting(members)) {
        // Collect a finished joint prediction.
        if (m_jointJob.valid() && m_jointJob.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            std::vector<trajectory::Path> paths = m_jointJob.get();
            for (size_t k = 0; k < m_jointSlots.size() && k < paths.size(); ++k) {
                if (m_jointSlots[k] >= 0) {
                    m_predictions[m_jointSlots[k]] = std::move(paths[k]);
                } else if (m_jointPlannerKey == m_ui.plannerKey) {
                    m_ui.plannerJointPath = std::move(paths[k]);
                    m_ui.plannerJointValid = true;
                }
            }
        }
        if (!m_jointJob.valid() && (m_settings.showPredictions || m_ui.plannerOpen)) {
            m_jointSlots = slots;
            double duration = m_settings.predictionTime;
            if (m_ui.plannerOpen && m_ui.plannerValid) {
                gravity::Member planned{m_ui.plannerConstants, m_ui.plannerState, m_ui.launchMass, 0.1 * m_ui.launchRadius};
                planned.s.t = m_bodies.time;
                members.push_back(planned);
                m_jointSlots.push_back(-1);
                m_jointPlannerKey = m_ui.plannerKey;
                duration = std::max(duration, 1500.0);
            }
            m_jointJob = std::async(std::launch::async, [spin, members, duration]() {
                return gravity::predict(spin, members, duration, 0.3, 3000);
            });
        }
        if (m_ui.plannerKey != m_jointPlannerKey) m_ui.plannerJointValid = false;
        return;
    }

    // Without interactions each body follows its own geodesic: a couple of bodies per
    // frame, so long predictions never stall a frame.
    m_ui.plannerJointValid = false;
    if (bodies.empty() || !m_settings.showPredictions) return;
    for (int k = 0; k < 2; ++k) {
        m_predictNext = (m_predictNext + 1) % static_cast<int>(bodies.size());
        const Body& b = bodies[m_predictNext];
        m_predictions[b.slot] = trajectory::predict(spin, b.constants, b.state, m_settings.predictionTime, 0.3, 3000);
    }
}

void App::buildPaths() {
    std::vector<PathLine> lines;
    auto toVec = [](const trajectory::Point& p) {
        return glm::vec3(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
    };
    const int riding = ridingIndex();
    if (m_settings.showPaths) {
        for (size_t i = 0; i < m_bodies.bodies.size(); ++i) {
            const Body& b = m_bodies.bodies[i];
            if (static_cast<int>(i) == riding) continue;  // you are on it
            const glm::vec3 color = displayColor(b) * 1.2f;
            PathLine trail;
            trail.color = color;
            trail.startOpacity = 0.03f;
            trail.endOpacity = 0.9f;
            const size_t stride = std::max<size_t>(1, b.trail.size() / 400);
            for (size_t j = 0; j < b.trail.size(); j += stride) trail.points.push_back(b.trail[j]);
            auto p = kerr::toCartesian(m_settings.spin, b.state.r, b.state.theta, b.state.phi);
            trail.points.push_back(toVec(p));
            lines.push_back(std::move(trail));

            const trajectory::Path& future = m_predictions[b.slot];
            if (m_settings.showPredictions && future.points.size() > 1) {
                PathLine predicted;
                predicted.color = color;
                predicted.dashed = true;
                predicted.startOpacity = 0.8f;
                predicted.endOpacity = 0.1f;
                const size_t step = std::max<size_t>(1, future.points.size() / 400);
                for (size_t j = 0; j < future.points.size(); j += step) predicted.points.push_back(toVec(future.points[j]));
                lines.push_back(std::move(predicted));
            }
        }
    }
    if (m_ui.plannerOpen && m_ui.plannerValid) {
        PathLine planned;
        planned.color = glm::vec3(0.3f, 0.85f, 1.0f);
        planned.dashed = true;
        planned.startOpacity = 1.0f;
        planned.endOpacity = 0.25f;
        const trajectory::Path& path = m_ui.plannerJointValid ? m_ui.plannerJointPath : m_ui.plannerPath;
        const size_t step = std::max<size_t>(1, path.points.size() / 600);
        for (size_t j = 0; j < path.points.size(); j += step) planned.points.push_back(toVec(path.points[j]));
        lines.push_back(std::move(planned));
    }
    m_renderer.setPaths(lines, m_settings.spin);
}

void App::inspectPixel(double x, double y) {
    int ww = 1, wh = 1;
    glfwGetWindowSize(m_window, &ww, &wh);
    const float t = std::tan(0.5f * glm::radians(m_camera.fovDeg));
    const float ndcX = static_cast<float>((2.0 * x - ww) / wh), ndcY = static_cast<float>((wh - 2.0 * y) / wh);
    const glm::vec3 dir = glm::normalize(m_camera.basis() * glm::vec3(ndcX * t, ndcY * t, 1.0f));
    const double look[3] = {dir.x, dir.y, dir.z};

    const double spin = m_settings.spin;
    trajectory::RayScene scene;
    scene.disk = m_settings.diskEnabled;
    scene.diskPrograde = m_settings.diskPrograde;
    scene.diskInner = m_settings.diskInnerAtIsco ? kerr::iscoRadius(spin, m_settings.diskPrograde)
                                                 : std::max<double>(m_settings.diskInner, kerr::horizonRadius(spin) + 0.05);
    scene.diskOuter = std::max<double>(m_settings.diskOuter, scene.diskInner + 0.5);
    scene.escapeRadius = std::max<double>(m_settings.escapeRadius, 1.5 * m_camera.r);
    std::vector<int> owner;
    const int riding = ridingIndex();
    for (size_t i = 0; i < m_bodies.bodies.size(); ++i) {
        if (static_cast<int>(i) == riding) continue;
        const Body& b = m_bodies.bodies[i];
        scene.spheres.push_back({kerr::toCartesian(spin, b.state.r, b.state.theta, b.state.phi), b.radius});
        owner.push_back(static_cast<int>(i));
    }

    m_ui.inspect = trajectory::traceRay(spin, m_camera.r, m_camera.theta, m_camera.phi, m_camera.frame(spin), look,
                                        scene);
    m_ui.inspectValid = true;
    m_ui.inspectorOpen = true;
    m_ui.inspectX = static_cast<float>(x);
    m_ui.inspectY = static_cast<float>(y);
    m_ui.inspectBody.clear();
    if (m_ui.inspect.body >= 0) m_ui.inspectBody = m_bodies.bodies[owner[m_ui.inspect.body]].name;
}

void App::applyScene(int scene) {
    Settings& s = m_settings;
    const Settings defaults;
    m_camera.reset();
    s.diskEnabled = true;
    s.diskPrograde = defaults.diskPrograde;
    s.diskInnerAtIsco = true;
    s.diskOuter = defaults.diskOuter;
    s.diskTemperature = defaults.diskTemperature;
    s.diskBrightness = defaults.diskBrightness;
    s.diskOpacity = defaults.diskOpacity;
    s.diskTurbulence = defaults.diskTurbulence;
    s.relativity = defaults.relativity;
    s.jetsEnabled = false;
    s.jetPower = defaults.jetPower;
    s.spin = defaults.spin;
    auto view = [&](double r, double inclinationDeg) {
        m_camera.r = r;
        m_camera.theta = inclinationDeg * kerr::kPi / 180.0;
    };
    switch (scene) {
        case SceneM87:
            // EHT: about 17 degrees from the spin axis; a powerful jet nearly toward us.
            s.spin = 0.94f;
            view(38.0, 17.0);
            s.diskTemperature = 5200.0f;
            s.diskBrightness = 5.0f;
            s.diskTurbulence = 0.85f;
            s.jetsEnabled = true;
            s.jetPower = 1.5f;
            break;
        case SceneSgrA:
            // Dim and turbulent; EHT images favour a fairly face-on view.
            s.spin = 0.9f;
            view(30.0, 35.0);
            s.diskTemperature = 3800.0f;
            s.diskBrightness = 3.0f;
            s.diskOpacity = 0.75f;
            s.diskTurbulence = 1.0f;
            break;
        case SceneGargantua:
            // Interstellar: near-maximal spin, nearly edge-on, frequency shifts left out.
            s.spin = 0.999f;
            view(26.0, 86.0);
            s.diskTemperature = 4600.0f;
            s.diskBrightness = 5.0f;
            s.diskOuter = 15.0f;
            s.diskTurbulence = 0.55f;
            s.relativity = 2;
            break;
        case SceneCygnusX1: {
            // About 27 degrees inclination, a hot inner disk, and a blue supergiant companion.
            s.spin = 0.998f;
            view(50.0, 27.0);
            s.diskTemperature = 16000.0f;
            s.diskBrightness = 3.0f;
            s.jetsEnabled = true;
            s.jetPower = 0.4f;
            if (m_bodies.addCircular(s.spin, 60.0, true, kerr::kPi * 0.8, BodyKind::Star, "Companion")) {
                m_bodies.bodies.back().temperature = 28000.0f;
                m_bodies.bodies.back().radius = 6.0f;
            }
            break;
        }
        default:
            break;
    }
    m_ui.cameraStatus = std::string("Scene: ") + kSceneNames[std::clamp(scene, 0, SceneCount - 1)];
}

void App::addEinsteinRing() {
    // Directly behind the hole as seen from the camera, at rest relative to the
    // local ZAMO: it falls almost straight in, staying lined up for a while.
    const double spin = m_settings.spin;
    const auto c = kerr::toCartesian(spin, m_camera.r, m_camera.theta, m_camera.phi);
    const double scale = -40.0 / std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    const auto [r, theta, phi] = kerr::fromCartesian(spin, c[0] * scale, c[1] * scale, c[2] * scale);
    const double v[3] = {0.0, 0.0, 0.0};
    if (m_bodies.addFromZamo(spin, r, theta, phi, v, BodyKind::Star, "Ring star")) {
        m_bodies.bodies.back().temperature = 9000.0f;
        m_bodies.bodies.back().radius = 1.6f;
        m_bodies.logEvent("A star sits exactly behind the hole: its light forms an Einstein ring.");
    }
}

void App::advanceSimulation(double dt) {
    const double spin = m_settings.spin;
    m_bodies.update(spin, dt, m_settings.bodyGravity);
    for (const auto& batch : m_bodies.pendingDebris) m_particles.spawn(batch.particles, m_bodies.time, batch.existedBefore);
    m_bodies.pendingDebris.clear();
    GasParams gas;
    gas.physics = m_settings.gasPhysics;
    gas.collisionRate = m_settings.gasCollisions;
    gas.coolingTime = m_settings.gasCoolingTime;
    gas.floorTemperature = m_settings.gasFloorTemperature;
    if (m_settings.bodyGravity) {
        for (const Body& b : m_bodies.bodies) {
            if (b.mass < 1e-5) continue;  // too light to stir the gas
            const auto p = kerr::toCartesian(spin, b.state.r, b.state.theta, b.state.phi);
            gas.bodies.emplace_back(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]),
                                    static_cast<float>(b.mass));
        }
    }
    // Light delay for the gas is measured from a camera outside the horizon.
    gas.lightDelay = m_settings.lightDelay && !m_camera.insideHorizon(spin);
    const auto cam = kerr::toCartesian(spin, m_camera.r, m_camera.theta, m_camera.phi);
    gas.cameraPosition = glm::vec3(static_cast<float>(cam[0]), static_cast<float>(cam[1]), static_cast<float>(cam[2]));
    gas.cameraR = static_cast<float>(m_camera.r);
    gas.time = m_bodies.time;
    m_particles.update(spin, dt, gas);
    const bool wasInside = m_camera.insideHorizon(spin), wasEnded = m_camera.ended;
    m_camera.advance(spin, dt);
    if (m_camera.insideHorizon(spin) && !wasInside && m_camera.motion != CameraMotion::Hover) {
        m_ui.cameraStatus = "You crossed the event horizon. Nothing marked the moment, but there is no way "
                            "back: from here every path, even light's, leads inward.";
    }
    if (m_camera.ended && !wasEnded) {
        m_ui.cameraStatus = spin > 0.01
            ? "The fall ends just above the inner horizon. There the Kerr interior stops being trustworthy "
              "(infalling light piles up without limit), so the simulation stops. Press R to start over."
            : "The fall ends just above the singularity, where tides crush everything. Press R to start over.";
        m_settings.paused = true;
    }
}

bool App::startRecording(const std::string& path, int width, int height, int fps, float quality, int cameraPath,
                         double seconds) {
    ensureParent(path);
    std::string error;
    if (!m_recording.encoder.open(path, width, height, fps, error)) {
        m_ui.videoStatus = "Could not start: " + error;
        std::fprintf(stderr, "Recording failed: %s\n", error.c_str());
        return false;
    }
    m_recording.active = true;
    m_recording.path = path;
    m_recording.width = width;
    m_recording.height = height;
    m_recording.fps = fps;
    m_recording.quality = quality;
    m_recording.cameraPath = cameraPath;
    m_recording.frame = 0;
    m_recording.total = std::max(1, static_cast<int>(std::lround(seconds * fps)));
    m_recording.startPhi = m_camera.phi;
    m_recording.startR = m_camera.r;
    glfwSwapInterval(0);  // render as fast as possible while recording
    m_ui.videoStatus = "Recording " + path;
    std::printf("Recording %d frames to %s\n", m_recording.total, path.c_str());
    return true;
}

void App::recordFrame() {
    int w = 0, h = 0;
    m_renderer.readOutput(m_recording.pixels, w, h);
    std::string error;
    if (!m_recording.encoder.write(m_recording.pixels.data(), error)) {
        m_ui.videoStatus = "Recording failed: " + error;
        std::fprintf(stderr, "%s\n", m_ui.videoStatus.c_str());
        stopRecording();
        return;
    }
    if (++m_recording.frame >= m_recording.total) stopRecording();
}

void App::stopRecording() {
    if (!m_recording.active) return;
    m_recording.active = false;
    std::string error;
    if (m_recording.encoder.close(error)) {
        char text[512];
        std::snprintf(text, sizeof(text), "Saved %s (%d frames, %.1f s)", m_recording.path.c_str(),
                      m_recording.frame, static_cast<double>(m_recording.frame) / m_recording.fps);
        m_ui.videoStatus = text;
    } else {
        m_ui.videoStatus = "Failed to finish the video: " + error;
    }
    std::printf("%s\n", m_ui.videoStatus.c_str());
    glfwSwapInterval(m_vsync ? 1 : 0);
}

int App::run(const LaunchOptions& options, const std::filesystem::path& shaderDir) {
    if (!init(options, shaderDir)) {
        shutdown();
        return 1;
    }
    const bool screenshotMode = !options.screenshotPath.empty();
    const bool recordMode = !options.recordPath.empty();
    const bool headless = screenshotMode || recordMode;
    if (recordMode && !startRecording(options.recordPath, options.width, options.height, options.recordFps,
                                      m_settings.renderScale, options.recordPathMode, options.recordSeconds)) {
        shutdown();
        return 1;
    }

    int frame = 0;
    double last = glfwGetTime();
    double lastReloadCheck = last;
    int exitCode = 0;

    while (!glfwWindowShouldClose(m_window)) {
        glfwPollEvents();
        int fbw = 0, fbh = 0;
        glfwGetFramebufferSize(m_window, &fbw, &fbh);
        if (headless) {
            fbw = options.width;
            fbh = options.height;
        }
        if (fbw == 0 || fbh == 0) {  // minimized
            glfwWaitEvents();
            continue;
        }

        const double now = glfwGetTime();
        // Recording and screenshots step time by fixed frames, so videos play smoothly.
        float dt = static_cast<float>(std::min(now - last, 0.1));
        if (m_recording.active) dt = 1.0f / static_cast<float>(m_recording.fps);
        else if (headless) dt = 1.0f / 60.0f;
        last = now;
        m_fps = m_fps <= 0.0f ? 1.0f / std::max(dt, 1e-4f) : 0.95f * m_fps + 0.05f / std::max(dt, 1e-4f);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        if (!headless) handleInput(dt);

        // Simulation
        if (m_settings.spin != m_lastSpin) {
            m_bodies.respin(m_lastSpin, m_settings.spin);
            m_particles.respin(m_lastSpin, m_settings.spin);
            m_camera.respin(m_lastSpin, m_settings.spin);
            m_lastSpin = m_settings.spin;
            m_camera.r = std::max(m_camera.r, kerr::horizonRadius(m_settings.spin) + 0.3);
        }
        const bool advance = !m_settings.paused && !screenshotMode;
        const double simDt = advance ? dt * m_settings.timeScale * approachSlowdown() : 0.0;
        advanceSimulation(simDt);
        syncRide(simDt);
        if (m_camera.lookAhead && m_camera.motion != CameraMotion::Hover) m_camera.faceVelocity(m_settings.spin);
        updatePredictions();

        if (m_recording.active && m_camera.motion == CameraMotion::Hover) {
            const double progress = static_cast<double>(m_recording.frame) / m_recording.total;
            if (m_recording.cameraPath == RecordOrbit) {
                m_camera.phi = m_recording.startPhi + 2.0 * kerr::kPi * progress;
            } else if (m_recording.cameraPath == RecordPushIn) {
                m_camera.r = std::max(m_recording.startR * std::pow(0.3, progress),
                                      kerr::horizonRadius(m_settings.spin) + 1.0);
            }
        } else if (m_settings.autoOrbit && !headless && m_camera.motion == CameraMotion::Hover) {
            m_camera.orbit(m_settings.autoOrbitSpeed * dt, 0.0);
        }
        if (m_settings.vsync != m_vsync && !headless && !m_recording.active) {
            m_vsync = m_settings.vsync;
            glfwSwapInterval(m_vsync ? 1 : 0);
        }
        if (now - lastReloadCheck > 0.5) {
            m_renderer.hotReload();
            m_particles.hotReload();
            lastReloadCheck = now;
        }

        // Render
        m_renderer.uploadHistory(m_bodies.history);
        buildPaths();
        GpuBodies gpuBodies;
        m_bodies.fillGpu(m_settings.spin, m_settings.lightDelay, m_camera.r, gpuBodies);
        FrameInputs in;
        in.settings = &m_settings;
        in.camera = &m_camera;
        in.bodies = &gpuBodies;
        in.particles = &m_particles;
        in.skipBody = ridingIndex();
        in.spin = m_settings.spin;
        in.time = m_bodies.time;
        in.width = m_recording.active ? m_recording.width : fbw;
        in.height = m_recording.active ? m_recording.height : fbh;
        if (m_recording.active) in.renderScale = m_recording.quality;
        m_renderer.render(in);
        if (m_recording.active) recordFrame();
        if (!headless) m_renderer.present(fbw, fbh);

        UiContext ui{m_settings, m_camera, m_bodies, m_renderer, m_particles, m_ui, m_fps, m_gpuName.c_str()};
        ui.recording = m_recording.active;
        ui.recordFrame = m_recording.frame;
        ui.recordTotal = m_recording.total;
        ui.predictions = &m_predictions;
        ui.slowdown = approachSlowdown();
        if (!headless || options.screenshotUi) drawUi(ui);
        ImGui::Render();
        if (headless) {
            // Draw the interface into the image the screenshot reads.
            m_renderer.bindOutput();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        } else {
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        }

        if (ui.screenshotRequested) saveScreenshot(timestampedPath("screenshots", ".png"));
        if (ui.cameraMotionRequested >= 0) setCameraMotion(static_cast<CameraMotion>(ui.cameraMotionRequested));
        if (ui.rideRequested >= 0) rideWith(ui.rideRequested);
        if (ui.sceneRequested >= 0) applyScene(ui.sceneRequested);
        if (ui.einsteinRingRequested) addEinsteinRing();
        if (ui.stopRecordingRequested) stopRecording();
        if (ui.recordRequested && !m_recording.active) {
            const VideoSize& size = kVideoSizes[m_ui.videoSize];
            startRecording(timestampedPath("videos", ".mp4"), size.width, size.height, kVideoFps[m_ui.videoFps],
                           m_ui.videoQuality, m_ui.videoPath, m_ui.videoSeconds);
        }

        glfwSwapBuffers(m_window);

        if (recordMode && !m_recording.active) break;  // video finished
        if (screenshotMode && ++frame >= options.frames) {
            glFinish();
            if (!m_renderer.shaderErrors().empty() || !m_particles.errors().empty()) {
                std::fprintf(stderr, "%s%s\n", m_renderer.shaderErrors().c_str(), m_particles.errors().c_str());
                exitCode = 2;
            }
            saveScreenshot(options.screenshotPath);
            for (const Body& b : m_bodies.bodies) {
                const auto x = kerr::toCartesian(m_settings.spin, b.state.r, b.state.theta, b.state.phi);
                std::printf("  %-16s %-12s r = %8.3f  theta = %6.3f  E = %.6f  L = %+.5f  Q = %.5f  "
                            "at (%.2f, %.2f, %.2f)\n",
                            b.name.c_str(), kindName(b.kind), b.state.r, b.state.theta, b.constants.E,
                            b.constants.L, b.constants.Q, x[0], x[1], x[2]);
            }
            if (m_particles.active()) std::printf("  debris: %s\n", m_particles.describe().c_str());
            if (m_camera.motion != CameraMotion::Hover) {
                glm::vec3 v = m_camera.velocity(m_settings.spin);
                std::printf("  camera: r = %.3f, speed %.3f c, clock rate %.3f, proper time %.2f M%s\n", m_camera.r,
                            glm::length(v), m_camera.clockRate(m_settings.spin), m_camera.properTime,
                            m_camera.ended ? " (the fall has ended)"
                                           : m_camera.insideHorizon(m_settings.spin) ? " (inside the horizon)" : "");
            }
            for (const std::string& e : m_bodies.events) std::printf("  event: %s\n", e.c_str());
            std::printf("GPU time per frame: %.2f ms at %dx%d traced\n", m_renderer.gpuMilliseconds(),
                        m_renderer.traceWidth(), m_renderer.traceHeight());
            break;
        }
    }
    shutdown();
    return exitCode;
}
