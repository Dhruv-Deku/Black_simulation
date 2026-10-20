#include "app.h"

#include "blackbody.h"
#include "central.h"
#include "file_dialog.h"
#include "keyframes.h"
#include "gravity.h"
#include "kerr.h"
#include "scene_file.h"
#include "surface.h"

#include <glad/gl.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <map>
#include <sstream>

namespace {

// The gas grid reaches this far (M) while a star's unbound debris is flying out.
constexpr float kFarDebrisRadius = 400.0f;
// A tidal disruption's outflow is drawn out to this radius (M), and the flare lights it
// this brightly at its launch radius when it peaks (with the panel's brightness at 1).
constexpr double kOutflowReach = 450.0;
constexpr float kOutflowLightLevel = 1.1f;
constexpr int kOutflowRows = 2048;
// How much faster time runs while a disrupted star's debris is far out (flarePace).
constexpr double kFastForward = 6.0;
// Starquakes on a young neutron star come every this many turns on average (Vela's
// glitches come every ~3 years: ~1e9 turns), swing its magnetic poles by up to this many
// radians at first (exaggerated: the crust moves by centimetres), and ring down with
// this quality factor (the core's field drains the crust's oscillations within seconds).
constexpr double kQuakeTurns = 18.0;
constexpr double kQuakeSwing = 0.05;
constexpr double kCrustQ = 60.0;
// A glitch's flash of X-rays: how bright, and its fade time (M).
constexpr double kQuakeFlash = 1.2, kQuakeFlashFade = 6000.0;
// Brightness of a beacon neutron star's needle beams.
constexpr float kBeaconBeamGlow = 2.4f;
// The timeline takes a checkpoint every this much scene time (times the scene's time
// unit), within these memory budgets (bytes; older checkpoints are thinned out to stay
// inside), and replays for at most this long (seconds) per frame while seeking.
constexpr double kCheckpointInterval = 25.0;
constexpr size_t kTimelineMemory = size_t(256) << 20;
constexpr size_t kTimelineGpuMemory = size_t(512) << 20;
constexpr double kSeekBudget = 0.03;

const std::map<std::string, Preset>& presetNames() {
    static const std::map<std::string, Preset> names = {
        {"isco", Preset::IscoOrbit},        {"circular", Preset::StableOrbit},  {"retrograde", Preset::RetrogradeOrbit},
        {"precessing", Preset::Precessing}, {"zoom-whirl", Preset::ZoomWhirl},  {"polar", Preset::PolarOrbit},
        {"plunge", Preset::Plunge},         {"flyby", Preset::Flyby},           {"cluster", Preset::Cluster},
        {"planet", Preset::RockyPlanet},    {"gas-giant", Preset::GasGiant},    {"pulsar", Preset::Pulsar},
        {"tde", Preset::TidalDisruption},   {"cloud", Preset::GasCloud},        {"ring", Preset::DebrisRing},
        {"tde-flare", Preset::TidalDisruptionFlare},
        {"grazing", Preset::GrazingStar},   {"binary", Preset::BinaryStar},     {"moon", Preset::PlanetWithMoon},
        {"mini-hole", Preset::MiniBlackHole}, {"white-dwarf", Preset::WhiteDwarfCompanion},
        // Real companions of stars and planets.
        {"the-moon", Preset::Moon},         {"space-station", Preset::SpaceStation},
        {"geostationary", Preset::GeostationarySatellite}, {"asteroid-impact", Preset::AsteroidImpact},
        {"galilean-moons", Preset::GalileanMoons}, {"comet-impact", Preset::CometImpact},
        {"inner-planets", Preset::InnerPlanets}, {"hot-jupiter", Preset::HotJupiter},
        {"sungrazer", Preset::SungrazingComet}, {"proxima-planets", Preset::ProximaPlanets},
        {"giant-planets", Preset::GiantInnerPlanets}, {"wd-planet", Preset::WhiteDwarfPlanet},
        {"wd-asteroid", Preset::WhiteDwarfAsteroid}, {"cv-companion", Preset::CataclysmicCompanion},
        {"clump", Preset::InfallingClump},
        {"penrose", Preset::Penrose},
    };
    return names;
}

// Command-line name of a catalogue entry: "Sirius B" -> "sirius-b", "Sgr A*" -> "sgr-a".
std::string centralKey(const char* label) {
    std::string key;
    for (const char* c = label; *c; ++c) {
        if (*c == ' ') key += '-';
        else if (*c != '*') key += static_cast<char>(std::tolower(static_cast<unsigned char>(*c)));
    }
    return key;
}

void APIENTRY glDebugCallback(GLenum, GLenum type, GLuint, GLenum severity, GLsizei, const GLchar* message,
                              const void*) {
    if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return;
    if (type == GL_DEBUG_TYPE_ERROR || severity == GL_DEBUG_SEVERITY_HIGH) {
        std::fprintf(stderr, "[GL] %s\n", message);
    }
}

void scrollCallback(GLFWwindow* window, double, double dy) {
    if (auto* inbox = static_cast<WindowInbox*>(glfwGetWindowUserPointer(window))) inbox->scroll += dy;
}

void dropCallback(GLFWwindow* window, int count, const char** paths) {
    if (auto* inbox = static_cast<WindowInbox*>(glfwGetWindowUserPointer(window))) {
        for (int i = 0; i < count; ++i) inbox->dropped.emplace_back(paths[i]);
    }
}

// Paths to and from UTF-8 (GLFW's and the scene files' encoding).
std::filesystem::path pathFromUtf8(const std::string& text) {
    return std::filesystem::path(std::u8string(reinterpret_cast<const char8_t*>(text.data()), text.size()));
}

std::string utf8Of(const std::filesystem::path& path) {
    const std::u8string text = path.u8string();
    return std::string(text.begin(), text.end());
}

json::Value encodeVec3(const glm::vec3& v) {
    json::Value a = json::Value::array();
    for (int i = 0; i < 3; ++i) a.push(json::Value::number(v[i]));
    return a;
}

bool decodeVec3(const json::Value* x, glm::vec3& v) {
    if (!x || !x->isArray() || x->size() != 3) return false;
    for (int i = 0; i < 3; ++i) {
        if (!(*x)[static_cast<size_t>(i)].get(v[i])) return false;
    }
    return true;
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
        "  --central NAME         the object at the center:\n"
        "                           earth jupiter brown-dwarf red-dwarf sun red-giant white-dwarf\n"
        "                           sirius-b neutron-star pulsar magnetar stellar-black-hole sgr-a m87\n"
        "  --mass MSUN            its mass in solar masses (sets the physical units)\n"
        "  --radius-km KM         its radius in km\n"
        "  --compactness R        its radius in M (R/M); below 9/4 it collapses\n"
        "  --temperature K        its surface temperature\n"
        "  --field G              its magnetic field in gauss (neutron stars, white dwarfs)\n"
        "  --physical-units       readouts in km, AU, seconds and hours too\n"
        "  --spin A               spin J/M^2 (a/M for a black hole, in [0, 0.999])\n"
        "  --naked                let a black hole spin past 1 (to 1.5): a naked singularity\n"
        "  --distance R           camera distance in M (around a star or planet, in units of R/6)\n"
        "  --inclination DEG      camera angle from the spin axis (90 = edge-on)\n"
        "  --azimuth DEG          camera azimuth\n"
        "  --fov DEG              vertical field of view\n"
        "  --camera hover|orbit|fall\n"
        "                         hold the camera still, put it in orbit, or drop it\n"
        "  --ride N               ride along with body number N (0 = first added)\n"
        "  --look-ahead           while moving, face the direction of travel\n"
        "  --scene default|m87|sgra|gargantua|cygnus|tde\n"
        "                         famous black holes; tde: a star torn apart, its flare and outflow\n"
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
        "  --no-flare             a torn-apart star's returning debris powers no flare or outflow\n"
        "  --starquake            a young neutron star's crust slips at the start and rings\n"
        "  --disk-temp K          peak disk temperature in kelvin\n"
        "  --disk-brightness B    disk brightness multiplier\n"
        "  --debris-glow B        debris brightness multiplier\n"
        "  --debris-thickness K   debris opacity (low = translucent glow)\n"
        "  --preset NAME          add objects (may repeat):\n"
        "                           orbits:  isco circular retrograde precessing zoom-whirl\n"
        "                                    polar plunge flyby cluster\n"
        "                           objects: planet gas-giant pulsar tde tde-flare grazing binary moon\n"
        "                                    cloud ring einstein-ring mini-hole white-dwarf penrose\n"
        "  --advance T            run the simulation for T M before starting\n"
        "  --scale S              ray tracing resolution scale (0.25 to 2)\n"
        "  --upscale bilinear|fsr how a smaller traced image is enlarged (default fsr: FSR 1)\n"
        "  --dynamic-resolution   the traced resolution follows a frame-time target\n"
        "  --no-disk              hide the accretion disk\n"
        "  --heatmap              show integration steps per pixel instead of the image\n"
        "  --debug-directions     show each ray's final direction as a color\n"
        "  --screenshot PATH      render without the UI, save a PNG and exit\n"
        "  --frames N             frames to render before that screenshot (default 4)\n"
        "  --ui                   include the interface in that screenshot\n"
        "  --record PATH          record an MP4 video (at --size) and exit\n"
        "  --duration S           video length in seconds (default 5)\n"
        "  --fps N                video frame rate (default 30)\n"
        "  --camera-path flown|orbit|push-in|keyframes\n"
        "                         how the camera moves during the video\n"
        "  --merger               two neutron stars merge (GW170817): inspiral, merger, black hole,\n"
        "                         gamma-ray burst, kilonova\n"
        "  --merger-view cinematic|simulation|split|film\n"
        "                         film-like view, simulation from above, or both side by side\n"
        "  --merger-particles N   SPH particles per star (default 12000)\n"
        "  --merger-masses M1,M2  the stars' masses in Suns (default 1.35,1.35)\n"
        "  --merger-collapse MS   the remnant collapses this long after the merger (default 15)\n"
        "  --merger-advance S     run the merger for S seconds of viewing before starting\n"
        "  --distance KM          around the merger, the camera's distance in km\n"
        "  --no-sound             no gravitational-wave sound (live or in videos)\n"
        "  --sound                the scene's sound: bodies, pulsars, the tide's rumble, waves\n"
        "  --no-captions          no captions in the merger's views\n"
        "  --scene-file PATH      open a saved scene (.bhscene) first; other options apply on top\n"
        "  --save-scene PATH      save the scene to PATH once set up (after --advance)\n"
        "  --seek T               after --advance, go back to time T on the timeline (replayed)\n"
        "  --no-timeline          do not record a timeline\n"
        "  --keyframe R,INCL,AZ[,S]\n"
        "                         a camera keyframe (may repeat): distance, inclination and azimuth\n"
        "                         as for --distance and the rest, reached S seconds after the last (4)\n"
        "  --eht m87|sgra|here    what the Event Horizon Telescope sees (here: this black hole)\n"
        "  --eht-inclination DEG  its spin axis to the line of sight (M87*: 163, Sgr A*: 30)\n"
        "  --eht-size N           its image's pixels across (default 160)\n"
        "  --keyframes-physical   the camera moves along the keyframes for real (aberration);\n"
        "                         --camera-path keyframes films along them\n");
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
            out.spin = static_cast<float>(std::clamp(v, 0.0, central::kMaxNakedSpin));
        } else if (arg == "--naked") {
            out.naked = true;
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
            int c = choice({"flown", "orbit", "push-in", "keyframes"});
            if (c < 0) return false;
            out.recordPathMode = c;
        } else if (arg == "--scene") {
            int c = choice({"default", "m87", "sgra", "gargantua", "cygnus", "tde"});
            if (c < 0) return false;
            out.scene = c;
        } else if (arg == "--central") {
            const char* name = value();
            if (!name) return false;
            const auto& list = central::catalog();
            int found = -1;
            for (size_t k = 0; k < list.size(); ++k) {
                if (centralKey(list[k].label) == name) found = static_cast<int>(k);
            }
            if (found < 0) {
                error = std::string("unknown central object: ") + name;
                return false;
            }
            out.central = found;
        } else if (arg == "--mass") {
            if (!number(v)) return false;
            out.massSolar = std::clamp(v, 1e-8, 1e11);
        } else if (arg == "--radius-km") {
            if (!number(v)) return false;
            out.radiusKm = std::max(v, 1e-3);
        } else if (arg == "--compactness") {
            if (!number(v)) return false;
            out.compactness = std::max(v, 0.5);
        } else if (arg == "--temperature") {
            if (!number(v)) return false;
            out.surfaceTemperature = std::clamp(v, 10.0, 1e8);
        } else if (arg == "--field") {
            if (!number(v)) return false;
            out.magneticField = std::clamp(v, 0.0, 1e16);
        } else if (arg == "--physical-units") {
            out.physicalUnits = true;
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
        } else if (arg == "--starquake") {
            out.starquake = true;
        } else if (arg == "--no-flare") {
            out.noFlare = true;
        } else if (arg == "--no-gravity") {
            out.noBodyGravity = true;
        } else if (arg == "--heatmap") {
            out.debugView = 1;
        } else if (arg == "--debug-directions") {
            out.debugView = 2;
        } else if (arg == "--ui") {
            out.screenshotUi = true;
        } else if (arg == "--merger") {
            out.merger = true;
        } else if (arg == "--merger-view") {
            int c = choice({"cinematic", "simulation", "split", "film"});
            if (c < 0) return false;
            out.mergerView = c;
            out.merger = true;
        } else if (arg == "--merger-particles") {
            if (!number(v)) return false;
            out.mergerParticles = std::clamp(static_cast<int>(v), 1000, 60000);
            out.merger = true;
        } else if (arg == "--merger-masses") {
            const char* s = value();
            double a = 0.0, b = 0.0;
            if (!s || std::sscanf(s, "%lf,%lf", &a, &b) != 2 || a < 1.0 || b < 1.0 || a > 2.0 || b > 2.0) {
                if (error.empty()) error = "--merger-masses expects M1,M2 in Suns (1 to 2), e.g. 1.4,1.3";
                return false;
            }
            out.mergerMass1 = a;
            out.mergerMass2 = b;
            out.merger = true;
        } else if (arg == "--merger-collapse") {
            if (!number(v)) return false;
            out.mergerCollapseMs = std::clamp(v, 1.0, 100.0);
            out.merger = true;
        } else if (arg == "--merger-advance") {
            if (!number(v)) return false;
            out.mergerAdvance = std::max(v, 0.0);
            out.merger = true;
        } else if (arg == "--no-sound") {
            out.noSound = true;
        } else if (arg == "--sound") {
            out.sound = true;
        } else if (arg == "--no-captions") {
            out.noCaptions = true;
        } else if (arg == "--scene-file") {
            const char* s = value();
            if (!s) return false;
            out.sceneFile = s;
        } else if (arg == "--save-scene") {
            const char* s = value();
            if (!s) return false;
            out.saveScene = s;
        } else if (arg == "--seek") {
            if (!number(v)) return false;
            out.seek = std::max(v, 0.0);
        } else if (arg == "--upscale") {
            int c = choice({"bilinear", "fsr"});
            if (c < 0) return false;
            out.upscaler = c;
        } else if (arg == "--dynamic-resolution") {
            out.dynamicResolution = true;
        } else if (arg == "--no-timeline") {
            out.noTimeline = true;
        } else if (arg == "--keyframe") {
            const char* s = value();
            LaunchOptions::Keyframe k{0.0, 0.0, 0.0, 4.0};
            if (!s || std::sscanf(s, "%lf,%lf,%lf,%lf", &k.distance, &k.inclinationDeg, &k.azimuthDeg, &k.seconds) < 3 ||
                k.distance <= 0.0) {
                if (error.empty()) error = "--keyframe expects R,INCL,AZ or R,INCL,AZ,SECONDS";
                return false;
            }
            out.keyframes.push_back(k);
        } else if (arg == "--eht") {
            int c = choice({"m87", "sgra", "here"});
            if (c < 0) return false;
            out.eht = c;
        } else if (arg == "--eht-inclination") {
            if (!number(v)) return false;
            out.ehtInclinationDeg = std::clamp(v, 1.0, 179.0);
        } else if (arg == "--eht-spin") {
            if (!number(v)) return false;
            out.ehtSpin = std::clamp(v, 0.0, 0.998);
        } else if (arg == "--eht-thickness") {
            if (!number(v)) return false;
            out.ehtHeight = std::clamp(v, 0.05, 2.0);
        } else if (arg == "--eht-temperature") {
            if (!number(v)) return false;
            out.ehtTemperature = std::clamp(v, 1e9, 1e14);
        } else if (arg == "--eht-size") {
            if (!number(v)) return false;
            out.ehtSize = std::clamp(static_cast<int>(v), 16, 1024);
        } else if (arg == "--keyframes-physical") {
            out.keyframesPhysical = true;
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
    glfwSetWindowUserPointer(m_window, &m_inbox);
    glfwSetScrollCallback(m_window, scrollCallback);
    glfwSetDropCallback(m_window, dropCallback);

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
    // The merger film's captions are a heavy italic, as the reference's are (Segoe UI Black
    // Italic, else Bold Italic).
    for (const char* captionFont : {"C:\\Windows\\Fonts\\seguibli.ttf", "C:\\Windows\\Fonts\\segoeuiz.ttf"}) {
        if (!std::filesystem::exists(captionFont)) continue;
        m_merger.setCaptionFont(io.Fonts->AddFontFromFileTTF(captionFont, 32.0f));
        break;
    }
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
    if (!m_merger.init(shaderDir, error)) {
        std::fprintf(stderr, "Merger shader errors:\n%s\n", error.c_str());
    }
    applyOptions(options);
    return true;
}

void App::applyOptions(const LaunchOptions& o) {
    // The timeline records in the window; for a picture or a video only if it is to seek.
    const bool headless = !o.screenshotPath.empty() || !o.recordPath.empty();
    m_timelineOn = !o.noTimeline && (!headless || o.seek.has_value());
    if (!o.sceneFile.empty()) loadScene(o.sceneFile);
    if (o.scene) applyScene(*o.scene);
    if (o.central) applyCentral(*o.central);
    central::Object& c = m_settings.central;
    if (o.massSolar) {
        // Same size in km, unless a white dwarf or neutron star follows its mass-radius relation.
        const double km = central::lengthKm(c.radius, c.massSolar);
        c.massSolar = *o.massSolar;
        double followed = 0.0;
        if (c.kind == central::Kind::WhiteDwarf) followed = central::whiteDwarfRadiusKm(c.massSolar);
        if (c.kind == central::Kind::NeutronStar) followed = central::neutronStarRadiusKm(c.massSolar);
        if (!c.hasHorizon()) c.radius = (followed > 0.0 ? followed : km) / (central::kKmPerSolarMass * c.massSolar);
    }
    if (o.radiusKm && !c.hasHorizon()) c.radius = *o.radiusKm / (central::kKmPerSolarMass * c.massSolar);
    if (o.compactness && !c.hasHorizon()) c.radius = *o.compactness;
    if (o.surfaceTemperature) c.surfaceTemperature = *o.surfaceTemperature;
    if (o.magneticField) c.magneticField = *o.magneticField;
    if (o.physicalUnits) m_settings.physicalUnits = true;
    if (o.naked && c.hasHorizon()) c.naked = true;
    if (o.spin) c.chi = c.naked ? *o.spin : std::min(*o.spin, 0.999f);
    syncCentral(0.0, 0.0f);
    if (o.distance) m_camera.r = *o.distance * sceneScale();
    if (o.inclinationDeg) m_camera.theta = std::clamp(*o.inclinationDeg, 1.0, 179.0) * kerr::kPi / 180.0;
    if (o.azimuthDeg) m_camera.phi = *o.azimuthDeg * kerr::kPi / 180.0;
    if (o.fov) m_camera.fovDeg = *o.fov;
    m_camera.r = std::max(m_camera.r, m_camera.minDistance(m_settings.spin));
    if (o.sky) m_settings.skyMode = *o.sky;
    if (o.renderScale) m_settings.renderScale = *o.renderScale;
    if (o.upscaler) m_settings.upscaler = *o.upscaler;
    if (o.dynamicResolution) m_settings.dynamicResolution = true;
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
    if (o.noFlare) m_settings.tdeFlare = false;
    if (o.jets) m_settings.jetsEnabled = true;
    if (o.sound) m_settings.sound = true;
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
    if (o.starquake) starquake();
    if (!o.keyframes.empty()) {
        const Camera saved = m_camera;
        m_track.keys.clear();
        for (const LaunchOptions::Keyframe& k : o.keyframes) {
            m_camera.motion = CameraMotion::Hover;
            m_camera.r = std::max(k.distance * sceneScale(), m_camera.minDistance(m_settings.spin));
            m_camera.theta = std::clamp(k.inclinationDeg, 1.0, 179.0) * kerr::kPi / 180.0;
            m_camera.phi = k.azimuthDeg * kerr::kPi / 180.0;
            m_camera.yaw = m_camera.pitch = m_camera.roll = 0.0f;
            if (!m_track.keys.empty()) m_track.keys.back().duration = std::max(k.seconds, 0.1);
            m_track.keys.push_back(keyframes::fromCamera(m_camera, m_settings.spin));
        }
        m_camera = saved;
        m_track.aim = true;  // each faces the center: so does the camera between them
        m_track.physical = o.keyframesPhysical;
        m_track.rebuild();
    }
    // --advance is in M around a black hole, in scene time (scale^1.5 M) around a bigger object.
    const double tScale = sceneScale() * std::sqrt(sceneScale());
    for (double t = 0.0; t < o.advance; t += 1.0) {
        const double dt = std::min(1.0, o.advance - t) * tScale;
        if (m_timelineOn) timelineFrame(0.0f, dt);
        else simulateFrame(dt, 0.0f);
    }
    if (o.seek && m_timelineOn && !m_timeline.empty()) {
        const size_t target = m_timeline.positionAt(*o.seek * tScale);
        while (!seekTowards(target, 1e9)) {
        }
        std::printf("Replayed the timeline back to t = %.6g M\n", m_bodies.time);
    }
    if (m_settings.flareCamera && !o.distance && m_camera.motion == CameraMotion::Hover) {
        m_camera.r = flareCameraDistance();
    }
    if (o.lookAhead) {
        m_camera.lookAhead = true;
        if (m_camera.motion != CameraMotion::Hover) m_camera.faceVelocity(m_settings.spin);
    }
    if (!o.saveScene.empty()) saveScene(o.saveScene);
    if (o.merger) {
        if (o.mergerParticles) m_merger.params.perStar = *o.mergerParticles;
        if (o.mergerCollapseMs) m_merger.params.collapseMs = *o.mergerCollapseMs;
        if (o.mergerMass1) m_merger.params.m1 = *o.mergerMass1;
        if (o.mergerMass2) m_merger.params.m2 = *o.mergerMass2;
        if (o.mergerView) m_merger.view = *o.mergerView;
        m_merger.sound = !o.noSound;
        m_merger.captions = !o.noCaptions;
        m_merger.start(m_camera);
        // Around the merger --distance is in km, and setting the view stops the camera's
        // following.
        if (o.distance) m_camera.r = *o.distance / merger::kKm;
        if (o.inclinationDeg) m_camera.theta = std::clamp(*o.inclinationDeg * kerr::kPi / 180.0, 0.02, kerr::kPi - 0.02);
        if (o.azimuthDeg) m_camera.phi = *o.azimuthDeg * kerr::kPi / 180.0;
        if (o.fov) m_camera.fovDeg = *o.fov;
        const double step = 1.0 / 30.0;
        if (m_merger.view == MergerScene::Film) m_merger.setFilmTime(o.mergerAdvance);
        else
            for (double t = 0.0; t < o.mergerAdvance; t += step) m_merger.update(step, false, m_camera, 1000, false);
        if (o.distance) m_merger.autoCamera = false;
    }
    if (o.eht >= 0) {
        m_ui.ehtOpen = true;
        m_ui.ehtTarget = o.eht;
        if (o.ehtInclinationDeg && o.eht < 2) m_ui.ehtInclination[o.eht] = static_cast<float>(*o.ehtInclinationDeg);
        if (o.ehtSize) m_ui.ehtResolution = *o.ehtSize <= 96 ? 0 : (*o.ehtSize <= 160 ? 1 : 2);
        if (o.ehtSize) m_ehtSizeOverride = *o.ehtSize;
        if (o.ehtHeight) m_ui.ehtHeight[o.eht] = static_cast<float>(*o.ehtHeight);
        if (o.ehtSpin) m_ui.ehtSpin = static_cast<float>(*o.ehtSpin);
        if (o.ehtTemperature) m_ui.ehtTemperature[o.eht] = static_cast<float>(*o.ehtTemperature);
    }
    if (o.inspect) {
        inspectPixel(o.inspect->first, o.inspect->second);
        const trajectory::RayReport& r = m_ui.inspect;
        std::printf("Inspected pixel (%.0f, %.0f): light came from %s\n", o.inspect->first, o.inspect->second,
                    trajectory::rayEndName(r.end));
        std::printf("  r_end = %.3f  g = %.3f  turns = %.2f  r_min = %.3f  delay = %.1f M  b = %.3f  q = %.3f\n",
                    r.rEnd, r.g, r.turns, r.rMin, r.delay, r.b, r.q);
        if (r.end == trajectory::RayEnd::Surface) {
            std::printf("  left the surface at theta = %.3f, %.1f degrees from the normal (its frame)\n", r.surfaceTheta,
                        std::acos(std::clamp(r.emissionCos, -1.0, 1.0)) * 180.0 / kerr::kPi);
        }
    }
}

void App::shutdown() {
    if (m_recording.active) stopRecording();
    m_sound.close();
    m_eht.shutdown();
    if (!m_timeline.empty()) {
        std::printf("Timeline: t = %.6g to %.6g M, %zu frames, %zu checkpoints (%zu after changes), %.1f MB + %.1f MB "
                    "on the GPU\n",
                    m_timeline.start(), m_timeline.end(), m_timeline.length(), m_timeline.checkpointCount(),
                    m_timeline.pinnedCount(), m_timeline.bytes() / 1048576.0, m_timeline.gpuBytes() / 1048576.0);
    }
    clearTimeline();
    m_merger.shutdown();
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
    const double minR = m_camera.minDistance(m_settings.spin);
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
    if (!io.WantCaptureMouse && m_inbox.scroll != 0.0 && !moving) {
        m_camera.zoom(std::pow(0.9, m_inbox.scroll), minR);
        m_settings.flareCamera = false;  // you set the distance now
    }
    m_inbox.scroll = 0.0;

    if (io.WantCaptureKeyboard) return;
    const bool shiftHeld = glfwGetKey(m_window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                           glfwGetKey(m_window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
    // Space plays or pauses (in the past too); Shift+Space plays the past backwards.
    if (pressed(GLFW_KEY_SPACE)) {
        if (shiftHeld && m_timelineOn) timelineCommand(2);
        else if (m_timelineOn && !m_timeline.empty()) timelineCommand(1);
        else m_settings.paused = !m_settings.paused;
    }
    if (pressed(GLFW_KEY_K)) addKeyframe();
    if (pressed(GLFW_KEY_T)) m_ui.ehtOpen = !m_ui.ehtOpen;
    if (pressed(GLFW_KEY_LEFT_BRACKET)) timelineCommand(5);
    if (pressed(GLFW_KEY_RIGHT_BRACKET)) timelineCommand(6);
    if (pressed(GLFW_KEY_TAB)) m_ui.showPanel = !m_ui.showPanel;
    if (pressed(GLFW_KEY_H)) m_ui.showHelp = !m_ui.showHelp;
    if (pressed(GLFW_KEY_R)) {
        m_camera.reset();
        m_settings.flareCamera = false;
    }
    if (pressed(GLFW_KEY_G)) m_settings.skyMode = m_settings.skyMode == 1 ? 0 : 1;
    if (pressed(GLFW_KEY_I)) m_settings.rayInspector = !m_settings.rayInspector;
    if (pressed(GLFW_KEY_F12)) saveScreenshot(timestampedPath("screenshots", ".png"));
    if (pressed(GLFW_KEY_ESCAPE) && m_recording.active) stopRecording();
    // Ctrl+S saves the scene (Ctrl+Shift+S: somewhere else), Ctrl+O opens one.
    const bool shift = glfwGetKey(m_window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                       glfwGetKey(m_window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
    const bool saveKey = pressed(GLFW_KEY_S), openKey = pressed(GLFW_KEY_O);
    if (ctrl && saveKey && !m_recording.active) saveSceneInteractive(shift);
    if (ctrl && openKey && !m_recording.active) openSceneInteractive();

    if (moving || ctrl) return;
    auto held = [&](int key) { return glfwGetKey(m_window, key) == GLFW_PRESS; };
    const double speed = 1.2 * dt;
    if (held(GLFW_KEY_W) || held(GLFW_KEY_S)) m_settings.flareCamera = false;
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
    const double x = m_camera.r / (4.0 * central::boundaryRadius(m_settings.spin, m_camera.surface));
    return std::clamp(x * x, 0.04, 1.0);
}

double App::flarePace() const {
    const Settings& s = m_settings;
    if (!s.flareFastForward || !s.tdeFlare || !s.central.hasHorizon() || m_bodies.tidalEvents.empty()) return 1.0;
    const tde::Event& e = m_bodies.tidalEvents.back();
    if (e.fallback.empty()) return 1.0;
    // In the time the camera sees (light from the center takes its distance to arrive):
    // from once the stream has left until just before it comes back.
    const double seen = m_bodies.time - (s.lightDelay ? m_camera.r : 0.0);
    const double from = e.time + 250.0, to = e.time + 0.9 * e.fallback.firstReturn;
    if (seen <= from || seen >= to) return 1.0;
    const double x = std::clamp(std::min(seen - from, to - seen) / 120.0, 0.0, 1.0);
    return 1.0 + (kFastForward - 1.0) * x * x * (3.0 - 2.0 * x);
}

double App::flareCameraDistance() const {
    const Settings& s = m_settings;
    if (m_bodies.tidalEvents.empty()) {
        // The star on its way in: keep it in view.
        for (const Body& b : m_bodies.bodies) {
            if (b.kind == BodyKind::TidalStar) return std::clamp(1.9 * b.state.r, 115.0, 160.0);
        }
        return m_camera.r;
    }
    const tde::Event& e = m_bodies.tidalEvents.back();
    const double seen = m_bodies.time - (s.lightDelay ? m_camera.r : 0.0);
    const double x = std::clamp((seen - e.time - 150.0) / 500.0, 0.0, 1.0);
    double target = 115.0 + 95.0 * x * x * (3.0 - 2.0 * x);
    if (s.tdeFlare && s.central.hasHorizon()) {
        const tde::Aftermath after(m_bodies.tidalEvents, outflowParams(s), s.central.massSolar, s.spin);
        if (after.active() && seen > after.start()) target = std::max(target, 1.15 * after.front(seen) + 95.0);
    }
    return std::min(target, 270.0);
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
    // Scene units: lengths scale with the scene, times with its dynamical time.
    const double scale = sceneScale(), tScale = scale * std::sqrt(scale);
    const double surface = m_bodies.surface;
    std::vector<gravity::Member> members;
    std::vector<int> slots;
    for (const Body& b : bodies) {
        members.push_back({b.constants, b.state, b.mass, 0.1 * b.radius});
        slots.push_back(b.slot);
    }
    // How far ahead, how finely and how far out. Around a star or planet, real orbits
    // take far longer than the scene's dynamical time: show at least one whole orbit,
    // its length from the energy (a = M / (2 (1 - E)), period 2 pi sqrt(a^3 / M)).
    struct Horizon {
        double duration, spacing, escape;
    };
    const bool realOrbits = !m_settings.central.hasHorizon();
    auto horizonFor = [&](const Body& b) {
        Horizon h{m_settings.predictionTime * tScale, 0.3 * scale, 2000.0 * scale};
        if (realOrbits && b.constants.E < 1.0) {
            const double a = 0.5 / (1.0 - b.constants.E);
            h.duration = std::max(h.duration, std::min(1.05 * 2.0 * kerr::kPi * a * std::sqrt(a), 3e5 * tScale));
            h.spacing = std::max(h.spacing, 2.0 * kerr::kPi * a / 1000.0);
            h.escape = std::max(h.escape, 2.5 * a);
        }
        return h;
    };

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
            Horizon all{m_settings.predictionTime * tScale, 0.3 * scale, 2000.0 * scale};
            // Kick-drift-kick steps: a dynamical time, or around a star or planet a
            // hundredth of the quickest orbit (a year needs no million steps).
            double maxStep = tScale, quickest = 1e300;
            for (const Body& b : bodies) {
                const Horizon h = horizonFor(b);
                all = {std::max(all.duration, h.duration), std::max(all.spacing, h.spacing),
                       std::max(all.escape, h.escape)};
                if (realOrbits && b.constants.E < 1.0) {
                    const double a = 0.5 / (1.0 - b.constants.E);
                    quickest = std::min(quickest, 2.0 * kerr::kPi * a * std::sqrt(a));
                }
            }
            if (quickest < 1e300) maxStep = std::max(maxStep, quickest / 100.0);
            double duration = all.duration;
            if (m_ui.plannerOpen && m_ui.plannerValid) {
                gravity::Member planned{m_ui.plannerConstants, m_ui.plannerState, m_ui.launchMass,
                                        0.1 * m_ui.launchRadius * scale};
                planned.s.t = m_bodies.time;
                members.push_back(planned);
                m_jointSlots.push_back(-1);
                m_jointPlannerKey = m_ui.plannerKey;
                duration = std::max(duration, 1500.0 * tScale);
            }
            m_jointJob = std::async(std::launch::async, [spin, members, duration, all, maxStep, surface]() {
                return gravity::predict(spin, members, duration, all.spacing, 3000, surface, all.escape, maxStep);
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
        const Horizon h = horizonFor(b);
        m_predictions[b.slot] = trajectory::predict(spin, b.constants, b.state, h.duration, h.spacing, 3000, surface,
                                                    h.escape);
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
    m_renderer.setPaths(lines, m_settings.spin, sceneScale(), !m_settings.central.hasHorizon());
}

void App::inspectPixel(double x, double y) {
    int ww = 1, wh = 1;
    glfwGetWindowSize(m_window, &ww, &wh);
    const float t = std::tan(0.5f * glm::radians(m_camera.fovDeg));
    const float ndcX = static_cast<float>((2.0 * x - ww) / wh), ndcY = static_cast<float>((wh - 2.0 * y) / wh);
    const glm::vec3 dir = glm::normalize(m_camera.basis() * glm::vec3(ndcX * t, ndcY * t, 1.0f));
    const double look[3] = {dir.x, dir.y, dir.z};

    const double spin = m_settings.spin;
    const CentralView view = centralView();
    trajectory::RayScene scene;
    scene.disk = m_settings.diskEnabled;
    scene.diskPrograde = m_settings.diskPrograde;
    scene.diskInner = view.diskInner;
    scene.diskOuter = view.diskOuter;
    scene.surface = view.surface;
    scene.surfaceOmega = view.omega;
    scene.escapeRadius = std::max<double>(m_settings.escapeRadius * sceneScale(), 1.5 * m_camera.r);
    std::vector<int> owner, lensOwner;
    const int riding = ridingIndex();
    for (size_t i = 0; i < m_bodies.bodies.size(); ++i) {
        if (static_cast<int>(i) == riding) continue;
        const Body& b = m_bodies.bodies[i];
        const auto pos = kerr::toCartesian(spin, b.state.r, b.state.theta, b.state.phi);
        if (b.kind == BodyKind::MiniBlackHole) {
            scene.lenses.push_back({pos, b.mass});
            lensOwner.push_back(static_cast<int>(i));
            continue;
        }
        scene.spheres.push_back({pos, b.radius});
        owner.push_back(static_cast<int>(i));
    }

    m_ui.inspect = trajectory::traceRay(spin, m_camera.r, m_camera.theta, m_camera.phi, m_camera.frame(spin), look,
                                        scene);
    m_ui.inspectValid = true;
    m_ui.inspectorOpen = true;
    m_ui.inspectX = static_cast<float>(x);
    m_ui.inspectY = static_cast<float>(y);
    m_ui.inspectBody.clear();
    if (m_ui.inspect.end == trajectory::RayEnd::Body && m_ui.inspect.body >= 0) {
        m_ui.inspectBody = m_bodies.bodies[owner[m_ui.inspect.body]].name;
    } else if (m_ui.inspect.end == trajectory::RayEnd::Lens && m_ui.inspect.body >= 0) {
        m_ui.inspectBody = m_bodies.bodies[lensOwner[m_ui.inspect.body]].name;
    } else if (m_ui.inspect.end == trajectory::RayEnd::Surface) {
        m_ui.inspectBody = m_settings.central.name;
    }
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
    // Undo what a star or neutron star left behind (applyCentral): its exposure, its
    // jets' speed, its pace.
    s.jetSpeed = defaults.jetSpeed;
    s.starBrightness = defaults.starBrightness;
    s.galaxyBrightness = defaults.galaxyBrightness;
    s.bloomStrength = defaults.bloomStrength;
    s.timeScale = defaults.timeScale;
    s.spin = defaults.spin;
    s.gasFloorTemperature = defaults.gasFloorTemperature;
    s.debrisGlowFloor = defaults.debrisGlowFloor;
    s.debrisFilaments = defaults.debrisFilaments;
    s.flareCamera = defaults.flareCamera;
    s.showPaths = defaults.showPaths;
    s.debrisBrightness = defaults.debrisBrightness;
    s.debrisOpacity = defaults.debrisOpacity;
    s.exposure = defaults.exposure;
    // Every famous scene has a black hole at the center.
    s.central = defaults.central;
    auto hole = [&](const char* name, double massSolar, float spin) {
        s.central.name = name;
        s.central.massSolar = massSolar;
        s.central.chi = spin;
        s.spin = spin;
    };
    auto view = [&](double r, double inclinationDeg) {
        m_camera.r = r;
        m_camera.theta = inclinationDeg * kerr::kPi / 180.0;
    };
    switch (scene) {
        case SceneM87:
            // EHT: about 17 degrees from the spin axis; a powerful jet nearly toward us.
            hole("M87*", 6.5e9, 0.94f);
            view(38.0, 17.0);
            s.diskTemperature = 5200.0f;
            s.diskBrightness = 5.0f;
            s.diskTurbulence = 0.85f;
            s.jetsEnabled = true;
            s.jetPower = 1.5f;
            break;
        case SceneSgrA:
            // Dim and turbulent; EHT images favour a fairly face-on view.
            hole("Sgr A*", 4.3e6, 0.9f);
            view(30.0, 35.0);
            s.diskTemperature = 3800.0f;
            s.diskBrightness = 3.0f;
            s.diskOpacity = 0.75f;
            s.diskTurbulence = 1.0f;
            break;
        case SceneGargantua:
            // Interstellar: near-maximal spin, nearly edge-on, frequency shifts left out.
            hole("Gargantua", 1e8, 0.999f);
            view(26.0, 86.0);
            s.diskTemperature = 4600.0f;
            s.diskBrightness = 5.0f;
            s.diskOuter = 15.0f;
            s.diskTurbulence = 0.55f;
            s.relativity = 2;
            break;
        case SceneCygnusX1: {
            // About 27 degrees inclination, a hot inner disk, and a blue supergiant companion.
            hole("Cygnus X-1", 21.2, 0.998f);
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
        case SceneTde: {
            // AT2019qiz: a hole of about a million Suns (Nicholl et al. 2020); its spin is
            // unknown. Quiet until it eats a star: no disk. A dark sky, as in ESO's film.
            hole("AT2019qiz", 1e6, 0.7f);
            view(150.0, 74.0);
            m_camera.fovDeg = 52.0f;
            s.diskEnabled = false;
            s.starBrightness = 0.2f;
            s.galaxyBrightness = 0.12f;
            s.showPaths = false;
            s.bloomStrength = 0.12f;
            s.timeScale = 40.0f;
            // The debris is cold (its pressure must not blow the stream apart) but glows
            // blue-white, lit by the star's leftover heat and then by the flare.
            s.gasFloorTemperature = 2500.0f;
            s.debrisGlowFloor = 11000.0f;
            s.debrisFilaments = true;
            s.flareCamera = true;
            s.debrisBrightness = 5.0f;
            s.debrisOpacity = 0.15f;
            m_bodies.clear();
            m_particles.clear();
            syncCentral(0.0, 0.0f);
            m_bodies.addPreset(s.spin, Preset::TidalDisruptionFlare);
            if (!m_bodies.bodies.empty()) {
                // Look across the plane of the orbit, square-on to the loop the returning
                // debris will make (the pericenter lies ~150 degrees on from the start).
                m_camera.phi = m_bodies.bodies.back().state.phi + 2.6 + 0.5 * kerr::kPi;
            }
            break;
        }
        default:
            break;
    }
    m_ui.cameraStatus = std::string("Scene: ") + kSceneNames[std::clamp(scene, 0, SceneCount - 1)];
    syncCentral(0.0, 0.0f);
}

void App::applyCentral(int index) {
    const auto& list = central::catalog();
    if (index < 0 || index >= static_cast<int>(list.size())) return;
    const central::CatalogEntry& e = list[static_cast<size_t>(index)];
    const Settings defaults;
    m_settings.central = e.object;
    m_settings.diskEnabled = e.disk;
    m_settings.diskInnerAtIsco = true;
    m_settings.diskOuter = defaults.diskOuter;
    m_settings.diskTemperature = defaults.diskTemperature;
    m_settings.diskBrightness = defaults.diskBrightness;
    // Neutron-star X-ray binaries launch jets at about half the speed of light (Sco X-1:
    // 0.45 c, Fomalont et al. 2001).
    m_settings.jetsEnabled = e.jets;
    m_settings.jetSpeed = e.jets && !e.object.hasHorizon() ? 0.5f : defaults.jetSpeed;
    m_settings.timeScale = e.timeScale > 0.0 ? static_cast<float>(e.timeScale) : defaults.timeScale;
    // Exposure: a camera set for a star's surface records no faint stars behind it (the
    // Sun's disk is ~1e10 times brighter per area than the starry sky), and its glare
    // blooms; around a sunlit planet the sky is dimmed less.
    const central::Kind kind = e.object.kind;
    const bool star = kind == central::Kind::Star || kind == central::Kind::WhiteDwarf;
    const bool planet = kind == central::Kind::RockyPlanet || kind == central::Kind::GasGiant;
    m_settings.starBrightness = star ? 0.04f : (planet ? 0.06f : defaults.starBrightness);
    m_settings.galaxyBrightness = star ? 0.15f : (planet ? 0.15f : defaults.galaxyBrightness);
    const bool neutron = kind == central::Kind::NeutronStar;
    m_settings.bloomStrength = star ? 0.1f : (neutron ? (e.object.beacon ? 0.06f : 0.14f) : defaults.bloomStrength);
    if (neutron) m_settings.starBrightness = e.object.beacon ? 0.3f : 0.6f;
    m_flashTime = -1e30;
    m_burstFuel = 0.0;
    m_quakes.clear();
    // Films show a spinning neutron star without its field's cage of lines.
    m_settings.fieldLines = !e.object.beacon;
    m_ui.compactnessOrigin = -1;
    // Rescale everything first, then put the camera where this object wants it.
    syncCentral(0.0, 0.0f);
    m_camera.reset();
    syncCentral(0.0, 0.0f);
    m_camera.r = std::max(e.cameraDistance * sceneScale(), m_camera.minDistance(m_settings.spin));
    // A white dwarf's disk spreads over tens of its radii (a third of the way to its
    // companion); the scene's default edge is set for black holes.
    if (e.diskOuter > 0.0 && !e.object.hasHorizon()) {
        m_settings.diskOuter = static_cast<float>(e.diskOuter * e.object.radius / sceneScale());
    }
    if (e.cameraInclinationDeg > 0.0) m_camera.theta = e.cameraInclinationDeg * kerr::kPi / 180.0;
    m_camera.roll = static_cast<float>(e.cameraRollDeg * kerr::kPi / 180.0);
    // A star's disk and its surface are both blackbodies: at the same temperature they
    // are equally bright. The black hole's disk brightness is a free, visual choice;
    // here it follows from the disk's temperature against the surface's (the surface is
    // drawn at kPhotosphereLevel = 0.42 of its color temperature's brightness).
    if (m_settings.diskPhysicalTemperature > 0.0) {
        const double surfaceT = central::displayTemperature(e.object.surfaceTemperature);
        m_settings.diskBrightness = static_cast<float>(0.42 * blackbody::luminance(m_settings.diskTemperature) /
                                                       blackbody::luminance(surfaceT));
    }
    m_ui.centralStatus = std::string(e.label) + ": " + e.description;
}

void App::triggerFlash(double strength, double fadeScene, const std::string& message) {
    const double scale = sceneScale();
    m_flashTime = m_bodies.time;
    m_flashStrength = strength;
    m_flashFade = fadeScene * scale * std::sqrt(scale);
    if (!message.empty()) {
        m_bodies.logEvent(message);
        m_ui.centralStatus = message;
    }
}

void App::starquake() {
    const central::Object& o = m_settings.central;
    if (o.kind != central::Kind::NeutronStar) return;
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    Quake q;
    q.time = m_bodies.time;
    q.l = uni(m_rng) < 0.6 ? 2 : 3;
    q.amplitude = static_cast<float>(kQuakeSwing * (0.6 + 0.6 * uni(m_rng)));
    // Any direction across the magnetic axis, in the star's rotating frame.
    const double tilt = o.magneticTiltDeg * kerr::kPi / 180.0;
    const glm::vec3 m(static_cast<float>(std::sin(tilt)), 0.0f, static_cast<float>(std::cos(tilt)));
    const glm::vec3 e1 = glm::normalize(glm::cross(m, std::abs(m.z) < 0.9f ? glm::vec3(0, 0, 1) : glm::vec3(1, 0, 0)));
    const glm::vec3 e2 = glm::cross(m, e1);
    const float angle = static_cast<float>(2.0 * kerr::kPi * uni(m_rng));
    q.direction = std::cos(angle) * e1 + std::sin(angle) * e2;
    m_quakes.push_back(q);
    if (m_quakes.size() > 4) m_quakes.erase(m_quakes.begin());
    char text[400];
    std::snprintf(text, sizeof(text),
                  "Starquake: the crust slipped (a glitch: the spin jumps by ~2 parts in a million) and rings at %.0f Hz "
                  "(l = %d), shaking the magnetic poles and the beams. Vela glitches every ~3 years; here every few "
                  "turns, to be seen.",
                  central::crustTorsionalFrequency(q.l), q.l);
    triggerFlash(kQuakeFlash, kQuakeFlashFade, text);
}

void App::syncCentral(double simDt, float realDt) {
    central::Object& o = m_settings.central;
    // Pressure gives out: a collapse (to a black hole, or a white dwarf to a neutron star).
    std::string message;
    if (central::collapseIfUnstable(o, message)) {
        m_ui.centralStatus = message;
        m_bodies.logEvent(message);
        m_screenFlash = 1.0f;
        m_ui.compactnessOrigin = -1;
        if (!o.hasHorizon()) triggerFlash(8.0, 30.0, "");  // a newborn neutron star blazes
    }
    // J is kept as the star is squeezed (only the spin slider stops at mass shedding),
    // so a collapse inherits it; the surface itself never turns faster than it can.
    m_settings.spin = static_cast<float>(o.metricSpin());

    const double scale = o.sceneScale();
    const double surface = o.hasHorizon() ? 0.0 : o.radius;
    m_bodies.surface = surface;
    m_bodies.surfaceOmega = o.hasHorizon() ? 0.0 : o.surfaceOmega();
    m_bodies.massSolar = m_settings.physicalUnits && o.hasHorizon() ? o.massSolar : 0.0;
    m_bodies.host = central::hostOf(o);
    if (o.hasHorizon()) m_settings.timeScale = std::min(m_settings.timeScale, 500.0f);  // the black hole's range
    m_bodies.centralMassSolar = o.massSolar;
    m_bodies.centralRadiusKm = o.hasHorizon() ? 0.0 : central::lengthKm(o.radius, o.massSolar);
    m_bodies.centralTemperature = o.surfaceTemperature;

    // A neutron star's or white dwarf's disk shines at the temperature its accretion
    // rate gives it (~1e7 K: X-rays; a few 1e4 K: ultraviolet), drawn at the compressed
    // color temperature used for hot surfaces. Black holes keep the chosen visual one.
    m_settings.diskPhysicalTemperature = 0.0;
    if ((o.kind == central::Kind::NeutronStar || o.kind == central::Kind::WhiteDwarf) && m_settings.diskEnabled) {
        const double inner = central::lengthKm(centralView().diskInner, o.massSolar);
        const double peak = central::diskPeakTemperature(
            o.massSolar, m_settings.accretionRate * central::eddingtonRate(o.massSolar), inner);
        m_settings.diskTemperature = static_cast<float>(central::displayTemperature(peak));
        m_settings.diskPhysicalTemperature = peak;
    }
    m_camera.surface = surface;
    m_camera.scale = scale;

    // Bodies, gas and camera keep their places and local velocities through a change of
    // spin; a change of scale only changes the units (and the camera stays put in them).
    if (m_settings.spin != m_lastSpin || scale != m_lastScale) {
        if (m_settings.spin != m_lastSpin) {
            m_bodies.respin(m_lastSpin, m_settings.spin);
            m_camera.respin(m_lastSpin, m_settings.spin);
        }
        m_particles.reframe(m_lastSpin, m_settings.spin, m_lastScale, scale, surface);
        if (scale != m_lastScale) {
            if (m_camera.motion == CameraMotion::Hover) m_camera.r *= scale / m_lastScale;
            m_bodies.setScale(scale, m_settings.spin);
        }
        m_lastSpin = m_settings.spin;
        m_lastScale = scale;
    }
    if (m_camera.motion == CameraMotion::Hover) m_camera.r = std::max(m_camera.r, m_camera.minDistance(m_settings.spin));

    // Thermonuclear flashes: hydrogen and helium piling up on a neutron star ignite in a
    // Type I X-ray burst (about 5 MeV per nucleon, a few percent of what the infall
    // released); on a white dwarf, a nova. Fuel arrives with the disk's accretion and
    // with anything that strikes the surface.
    const bool neutron = o.kind == central::Kind::NeutronStar, dwarf = o.kind == central::Kind::WhiteDwarf;
    const double tScale = scale * std::sqrt(scale);
    if (neutron || dwarf) {
        double fuel = 0.0;
        if (m_settings.diskEnabled && simDt > 0.0) {
            fuel += m_settings.accretionRate / 0.01 * simDt / ((neutron ? 600.0 : 3000.0) * tScale);
        }
        const double gas = m_particles.accreted();
        const double struck = static_cast<double>(m_bodies.impactCount);
        fuel += 0.4 * std::max(0.0, struck + 2.0 * gas - m_lastAccreted);
        m_lastAccreted = struck + 2.0 * gas;
        m_burstFuel += fuel;
        if (m_burstFuel >= 1.0) {
            m_burstFuel = 0.0;
            if (neutron) {
                triggerFlash(4.0, 30.0, "Type I X-ray burst: the helium layer on the surface ignited and burned across the "
                                        "star in seconds (~5 MeV per nucleon).");
            } else {
                triggerFlash(10.0, 300.0, "Nova: hydrogen accreted onto the white dwarf ignited in a thermonuclear runaway.");
            }
        }
    }
    if (o.magnetar && simDt > 0.0 &&
        std::uniform_real_distribution<double>(0.0, 1.0)(m_rng) < simDt / (1500.0 * tScale)) {
        triggerFlash(6.0, 20.0, "Magnetar flare: the crust cracked under the magnetic stress.");
    }
    // A young neutron star's crust is strained as its spin slows; now and then it gives
    // (a glitch) and rings in its torsional modes.
    if (neutron && o.beacon && m_settings.starquakes && simDt > 0.0) {
        const double period = 2.0 * kerr::kPi / std::max(std::abs(o.surfaceOmega()), 1e-12);
        if (std::uniform_real_distribution<double>(0.0, 1.0)(m_rng) < simDt / (kQuakeTurns * period)) starquake();
    }

    // White-light flares on active stars (red dwarfs like Proxima flare every few hours):
    // magnetic loops over a spot group reconnect and dump their energy into the
    // chromosphere, which flashes at ~9,500 K for minutes. At the star's real rate, in
    // its own time.
    const surface::Appearance look = surface::appearance(o);
    const double secondsPerM = central::timeSeconds(1.0, o.massSolar);
    if (look.flareInterval > 0.0 && m_settings.surfaceDetail && simDt > 0.0) {
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        if (uni(m_rng) < simDt * secondsPerM / look.flareInterval) {
            const surface::SpotPattern pattern = surface::makeSpots(
                look.spotStyle, look.activity, static_cast<unsigned>(std::hash<std::string>{}(o.name)));
            glm::vec3 where(0.0f, 0.0f, 1.0f);
            if (!pattern.regions.empty()) {
                const surface::ActiveRegion& region =
                    pattern.regions[std::min(pattern.regions.size() - 1, static_cast<size_t>(uni(m_rng) * pattern.regions.size()))];
                const glm::vec3 jitter(static_cast<float>(uni(m_rng) - 0.5), static_cast<float>(uni(m_rng) - 0.5),
                                       static_cast<float>(uni(m_rng) - 0.5));
                where = glm::normalize(region.direction + jitter * region.radius);
            }
            m_flares.push_back({where, m_bodies.time, static_cast<float>(0.5 + 2.5 * uni(m_rng) * uni(m_rng))});
            if (m_flares.size() > 4) m_flares.erase(m_flares.begin());
            m_bodies.logEvent("Stellar flare: magnetic loops over a spot group snapped and reconnected; the surface "
                              "below flashed at ~9,500 K.");
        }
    }
    if (look.flareInterval <= 0.0) m_flares.clear();
    m_screenFlash *= std::exp(-std::max(realDt, 0.0f) / 0.35f);
}

CentralView App::centralView() const {
    const Settings& s = m_settings;
    const central::Object& o = s.central;
    const double spin = s.spin, scale = sceneScale();
    CentralView v;
    v.kind = static_cast<int>(o.kind);
    v.detail = s.surfaceDetail;
    v.brightness = s.surfaceBrightness;
    const double isco = kerr::iscoRadius(spin, s.diskPrograde);
    if (o.hasHorizon()) {
        v.diskInner = s.diskInnerAtIsco ? isco : std::max<double>(s.diskInner * scale, kerr::horizonRadius(spin) + 0.05);
    } else {
        v.surface = o.radius;
        v.omega = o.surfaceOmega();
        v.phase = std::fmod(v.omega * m_bodies.time, 2.0 * kerr::kPi);
        v.turns = std::fmod(v.omega * m_bodies.time / (2.0 * kerr::kPi), 4096.0);
        v.azimuthShift = kerr::azimuthShift(spin, o.radius);
        v.temperature = static_cast<float>(central::displayTemperature(o.surfaceTemperature));
        const double tilt = o.magneticTiltDeg * kerr::kPi / 180.0;
        v.magneticAxis = glm::vec3(static_cast<float>(std::sin(tilt)), 0.0f, static_cast<float>(std::cos(tilt)));
        v.pulsar = o.pulsar;
        v.magnetar = o.magnetar;
        const bool neutron = o.kind == central::Kind::NeutronStar, dwarf = o.kind == central::Kind::WhiteDwarf;
        // Polar caps on a magnetic white dwarf, where the accretion columns land (a
        // neutron star's are set below, once we know whether gas reaches them).
        if (dwarf && o.magneticField > 1e5 && s.diskEnabled) {
            v.hotSpotTemperature = static_cast<float>(central::displayTemperature(o.surfaceTemperature * 3.0));
        }
        // How the surface looks, from its physics (src/surface.h): the convection cells'
        // size and lifetime from the pressure scale height, starspots in the latitudes and
        // groups of a real star, the corona's real profile (brightened, since a real corona
        // is a millionth as bright as the disk), an Earth-like sky at true scale.
        const surface::Appearance look = surface::appearance(o);
        v.look = static_cast<int>(look.look);
        v.planetMaps = s.planetMaps;
        v.filmColor = s.filmColor && (o.kind == central::Kind::Star || o.kind == central::Kind::WhiteDwarf ||
                                      o.kind == central::Kind::NeutronStar);
        v.radiusKm = central::lengthKm(o.radius, o.massSolar);
        // Jupiter's Great Red Spot starts on the day side (longitude 2 seed - 2.4).
        v.planetSeed = o.kind == central::Kind::GasGiant ? 0.85f : 3.7f;
        v.planetLight = 1.0f;
        if (look.look == surface::Look::BrownDwarf) v.planetSeed = 2.3f;
        // Planets and moons among the bodies are lit by the star at the center (its real
        // flux, relative to the Sun's at the Earth), or around a planet by the same Sun
        // that lights it.
        if (o.kind == central::Kind::RockyPlanet || o.kind == central::Kind::GasGiant) {
            v.bodyLight = 2;
            v.bodyLightScale = v.planetLight;
        } else {
            v.bodyLight = 1;
            v.bodyLightTemperature = static_cast<float>(central::displayTemperature(o.surfaceTemperature));
            v.bodyLightScale = central::fluxRelativeToEarth(v.radiusKm, o.surfaceTemperature,
                                                            central::lengthKm(1.0, o.massSolar));
        }
        if (s.corona) {
            v.atmosphere = look.atmosphere;
            if (look.corona) {
                v.coronaScale = static_cast<float>(surface::coronaScaleReal() * s.coronaBoost);
                v.chromosphereScale = 1.5f;
            }
            // In film colors a white dwarf glows inside blue wisps (its corona's shape).
            if (o.kind == central::Kind::WhiteDwarf && v.filmColor) {
                v.coronaScale = static_cast<float>(surface::coronaScaleReal() * s.coronaBoost * 4.0);
            }
        }
        // A giant's extended atmosphere (molecular layers, shed dust), in film colors.
        v.sdoMap = v.filmColor && s.surfaceDetail && look.look == surface::Look::Photosphere && look.granuleLane <= 0.3;
        if (v.filmColor && s.corona && look.look == surface::Look::Photosphere && look.granuleLane > 0.3) {
            v.envelope = 0.7f;
        }
        if (look.look == surface::Look::Photosphere) {
            const double secondsPerM = central::timeSeconds(1.0, o.massSolar);
            // Flares rise in about a minute and fade over ~10 (in the star's time).
            v.flareRise = 60.0 / secondsPerM;
            v.flareDecay = 600.0 / secondsPerM;
            for (const StarFlare& flare : m_flares) {
                if (m_bodies.time - flare.time > 12.0 * v.flareDecay) continue;
                v.flares.push_back(glm::vec4(flare.direction, flare.strength));
                v.flareAges.push_back(m_bodies.time - flare.time);
            }
            v.granuleFrequency = static_cast<float>(look.granuleFrequency);
            v.granuleContrast = static_cast<float>(look.granuleContrast);
            v.granuleLane = static_cast<float>(look.granuleLane);
            v.granuleRate = secondsPerM / look.granuleLifetime;
            v.granulePhase = static_cast<float>(std::fmod(m_bodies.time * v.granuleRate, 4096.0));
            v.spotDrop = static_cast<float>(look.spotDrop);
            // The same star always has the same spots.
            const unsigned seed = static_cast<unsigned>(std::hash<std::string>{}(o.name));
            const surface::SpotPattern pattern = surface::makeSpots(look.spotStyle, look.activity, seed);
            for (const surface::Spot& spot : pattern.spots) {
                v.spots.push_back(glm::vec4(spot.direction, spot.radius));
                v.spotUmbra.push_back(spot.umbra);
            }
            for (const surface::ActiveRegion& region : pattern.regions) {
                v.regions.push_back(glm::vec4(region.direction, region.radius));
            }
        }
        // Polar caps: where the field lines that carry the heat or the gas land. A
        // pulsar's open field lines are those that reach the light cylinder c / Omega,
        // landing within sin^2(theta) = R Omega of the pole; gas from a disk stopped at
        // r_m lands within sin^2(theta) = R / r_m. A magnetar's twisted field heats far
        // larger areas.
        if (neutron || dwarf) {
            double sin2 = o.radius * std::abs(v.omega);
            if (o.magnetar) sin2 = std::max(sin2, 0.06);
            sin2 = std::clamp(sin2, 0.002, 0.5);
            v.capCos = static_cast<float>(std::sqrt(1.0 - sin2));
        }
        // Magnetospheric truncation: the field stops the disk near half the Alfven radius,
        // and the gas rides the field lines down to the poles.
        double inner = s.diskInnerAtIsco ? std::max(isco, o.radius) : std::max<double>(s.diskInner * scale, o.radius);
        if ((neutron || dwarf) && s.magnetosphere && s.diskEnabled && o.magneticField > 0.0) {
            const double rKm = central::lengthKm(o.radius, o.massSolar);
            const double rA = central::alfvenRadiusKm(o.magneticField, rKm, o.massSolar,
                                                      s.accretionRate * central::eddingtonRate(o.massSolar)) /
                              (central::kKmPerSolarMass * o.massSolar);
            if (0.5 * rA > o.radius * 1.05) {
                v.funnelRadius = 0.5 * rA;
                inner = std::max(inner, v.funnelRadius);
            }
        }
        v.diskInner = inner;
        if (v.funnelRadius > 0.0) {
            const double sin2 = std::clamp(o.radius / v.funnelRadius, 0.002, 0.5);
            v.capCos = static_cast<float>(std::sqrt(1.0 - sin2));
        }
        // A neutron star's magnetic poles are hot where something heats them: gas
        // falling down the accretion columns (~1e7 K, X-ray pulsars), or a pulsar's
        // return currents striking its polar caps (a few times the surface's temperature).
        // An old neutron star fed at its equator, or a magnetar (its heat is in one
        // current-heated spot, drawn by the shader), has none.
        if (neutron && !o.magnetar) {
            if (v.funnelRadius > 0.0) v.hotSpotTemperature = static_cast<float>(central::displayTemperature(1e7));
            else if (o.pulsar) v.hotSpotTemperature = static_cast<float>(central::displayTemperature(2.5 * o.surfaceTemperature));
        }
        // The magnetosphere: closed field lines out to the light cylinder c / Omega, or to
        // where the disk takes over (none when the disk reaches the surface: an old
        // neutron star's weak field is buried under the gas), and drawn no farther than
        // where they fade. A magnetar's field is twisted by ~1 radian (Thompson,
        // Lyutikov and Kulkarni 2002) and threaded by a corona of pairs. Rotation-powered
        // pulsars, and radio-loud magnetars, beam from their poles; while gas pours onto
        // the poles (an accreting X-ray pulsar) the radio beams are quenched.
        if (neutron && s.surfaceDetail && o.magneticField > 0.0) {
            const double lightCylinder = 1.0 / std::max(std::abs(v.omega), 1e-12);
            constexpr double kFieldDrawnRadii = 7.0;
            double reach = std::min(lightCylinder, kFieldDrawnRadii * o.radius);
            if (v.funnelRadius > 0.0) reach = std::min(reach, v.funnelRadius);
            else if (s.diskEnabled) reach = std::min(reach, inner);
            const float glow = s.magnetosphereGlow;
            if (s.fieldLines && reach > o.radius) {
                v.fieldReach = reach;
                v.fieldGlow = (o.magnetar ? 2.5f : 1.6f) * glow;
                v.fieldTwist = o.magnetar ? 1.0f : 0.0f;
                v.pairGlow = o.magnetar ? 0.6f * glow : 0.0f;
            }
            if (v.filmColor && s.diskEnabled) v.flowGlow = 1.0f;
            // Gas pouring onto the poles usually quenches a pulsar's radio beams, but
            // transitional millisecond pulsars (PSR J1023+0038) pulse even while they accrete.
            if (s.pulsarBeams && (o.pulsar || o.magnetar)) {
                const double period = 2.0 * kerr::kPi / std::max(std::abs(v.omega), 1e-12) *
                                      central::timeSeconds(1.0, o.massSolar);
                v.beamWidth = static_cast<float>(central::beamHalfAngle(period));
                // The same power in a narrower beam is brighter (flux ~ 1 / its solid angle).
                constexpr double kReferenceBeam = 0.25;  // radians
                const double narrowing = std::clamp(std::pow(kReferenceBeam / v.beamWidth, 2.0), 0.15, 6.0);
                v.beamGlow = static_cast<float>(3.0 * narrowing) * glow;
                if (o.beacon) {
                    // Films show the bright needle at the beam's center: the core.
                    v.beamWidth = static_cast<float>(central::beamCoreHalfWidth(period));
                    v.beamGlow = kBeaconBeamGlow * glow;
                }
            }
            if (o.beacon) {
                // Its glow, and the starquakes ringing in its crust (each fading within a
                // few hundred oscillations, Q = kCrustQ).
                v.beacon = true;
                v.haloGlow = glow;
                const double secondsPerM = central::timeSeconds(1.0, o.massSolar);
                for (const Quake& q : m_quakes) {
                    const double f = central::crustTorsionalFrequency(q.l);
                    const double decay = kCrustQ / (kerr::kPi * f) / secondsPerM;
                    const double age = m_bodies.time - q.time;
                    if (age > 6.0 * decay) continue;
                    v.quakes.push_back(glm::vec4(static_cast<float>(age), q.amplitude,
                                                 static_cast<float>(2.0 * kerr::kPi * f * secondsPerM),
                                                 static_cast<float>(decay)));
                    v.quakeAxes.push_back(glm::vec4(q.direction, static_cast<float>(q.l)));
                }
            }
        }
        // Bodies that struck the surface, where they hit in the rotating frame (long faded
        // ones are left out).
        v.impactFade = 20.0 * scale * std::sqrt(scale);
        for (const Impact& impact : m_bodies.impacts) {
            const double age = m_bodies.time - impact.time;
            if (age > 20.0 * v.impactFade) continue;
            const double lon = impact.phi - v.azimuthShift - v.omega * impact.time;
            const double st = std::sin(impact.theta);
            const float strength = static_cast<float>(3.0 * std::clamp(std::sqrt(impact.energy / 0.05), 0.3, 4.0));
            v.impacts.push_back(glm::vec4(static_cast<float>(st * std::cos(lon)), static_cast<float>(st * std::sin(lon)),
                                          static_cast<float>(std::cos(impact.theta)), strength));
            v.impactAges.push_back(age);
        }
        if (v.impacts.size() > 8) {
            v.impacts.erase(v.impacts.begin(), v.impacts.end() - 8);
            v.impactAges.erase(v.impactAges.begin(), v.impactAges.end() - 8);
        }
    }
    v.impactFade = 20.0 * scale * std::sqrt(scale);
    v.diskOuter = std::max(s.diskOuter * scale, v.diskInner * 1.2);
    // A white dwarf fed by a disk is a cataclysmic variable: its red-dwarf companion
    // (0.15 solar masses, 0.19 solar radii) fills its Roche lobe, a ~2 hour orbit away,
    // and its gas falls from L1 to the disk along the ballistic path (computed once per
    // disk size).
    if (o.kind == central::Kind::WhiteDwarf && s.diskEnabled) {
        constexpr double kCompanionMass = 0.15, kCompanionRadius = 0.19;  // solar units
        const double q = kCompanionMass / o.massSolar;
        const double aKm = kCompanionRadius * central::kSunRadiusKm / central::rocheLobeFraction(q);
        const double a = aKm / (central::kKmPerSolarMass * o.massSolar);
        const double omega = std::sqrt((1.0 + q) / (a * a * a));
        // Starts beside the default view's white dwarf, a little toward the camera.
        constexpr double kStartPhase = 1.0;
        const double phase = std::fmod(kStartPhase + omega * m_bodies.time, 2.0 * kerr::kPi);
        const double rim = std::min(v.diskOuter / a, central::innerLagrangeFraction(q) * 0.95);
        if (q != m_streamQ || std::abs(rim - m_streamRim) > 1e-4) {
            double closest = 0.0;
            m_streamPath = central::rocheStream(q, rim, 16, closest);
            m_streamQ = q;
            m_streamRim = rim;
        }
        const double c = std::cos(phase), sn = std::sin(phase);
        v.binaryPhase = static_cast<float>(phase);
        if (!m_streamPath.empty()) {
            const auto& end = m_streamPath.back();
            v.spotPhase = static_cast<float>(phase + std::atan2(end[1], end[0]));
        }
        if (s.companionStar) {
            v.companion = glm::vec4(static_cast<float>(a * c), static_cast<float>(a * sn), 0.0f,
                                    static_cast<float>(central::rocheLobeFraction(q) * a));
            v.companionAxis = glm::vec3(static_cast<float>(c), static_cast<float>(sn), 0.0f);
            v.companionQ = static_cast<float>(q);
            for (const auto& p : m_streamPath) {
                v.stream.push_back(glm::vec3(static_cast<float>(a * (p[0] * c - p[1] * sn)),
                                             static_cast<float>(a * (p[0] * sn + p[1] * c)), 0.0f));
            }
            v.streamWidth = 0.0045 * a;
        }
    }
    v.flashAge = m_bodies.time - m_flashTime;
    v.flashStrength = m_flashStrength;
    v.flashFade = m_flashFade;
    return v;
}

FlareView App::flareView() const {
    FlareView fv;
    const Settings& s = m_settings;
    if (!s.tdeFlare || !s.central.hasHorizon() || m_bodies.tidalEvents.empty()) return fv;
    const tde::Aftermath after(m_bodies.tidalEvents, outflowParams(s), s.central.massSolar, s.spin);
    if (!after.active()) return fv;
    const tde::Event& ref = after.reference();
    const double now = m_bodies.time;
    // The camera sees the center as it was a light-travel time ago.
    const double seen = now - (s.lightDelay ? m_camera.r : 0.0);
    const double R0 = after.launchRadius();
    fv.launchRadius = R0;
    fv.speed = s.outflowSpeed;
    fv.axis = glm::vec3(static_cast<float>(ref.axis[0]), static_cast<float>(ref.axis[1]), static_cast<float>(ref.axis[2]));
    // The new disk, from the moment a couple of percent of the star has come back.
    const float lumSeen = static_cast<float>(after.luminosity(seen));
    const float presence = static_cast<float>(std::clamp(after.returned(seen) / 0.02, 0.0, 1.0));
    fv.diskGlow = s.flareDisk ? lumSeen * presence : 0.0f;
    fv.diskOuter = 1.3 * ref.circularization;
    fv.diskTemperature = 2600.0f + 1700.0f * std::pow(std::max(lumSeen, 0.0f), 0.25f);
    if (now <= after.start()) return fv;

    fv.active = true;
    fv.extent = std::min(after.front(now) * 1.15 + 10.0, kOutflowReach);
    fv.clock = now - ref.time;
    fv.start = after.start() - ref.time;
    fv.thickness = s.outflowThickness;
    fv.brightness = s.outflowBrightness * kOutflowLightLevel;
    fv.falloff = s.flareCinematic ? 0.7f : 2.0f;
    // TDE flares are hot (a few 1e4 K): white-blue, reddened a little by the wind's recession.
    fv.flareTemperature = 9500.0f;
    fv.glow = 0.25f;
    fv.envelope = 0.0015f;
    // The gas leaves the circularization radius with some of its orbital angular momentum
    // (here half the Keplerian sqrt(M R0)) and keeps it: it turns about the axis by
    // ell / (v R0) on its way out.
    fv.swirl = static_cast<float>(0.5 / (std::max(s.outflowSpeed, 0.01f) * std::sqrt(R0)));

    // The history, far enough back for the slowest gas drawn (the envelope, at a third of the
    // speed) and for light from the far side of the outflow still on its way to the camera.
    const double travel = tde::windTravelTime(1.0, 0.33 * s.outflowSpeed, R0, std::max(fv.extent, R0 + 1.0));
    const double window = std::clamp(std::max(travel, fv.clock - fv.start + 100.0) + m_camera.r + fv.extent + 100.0,
                                      500.0, 40000.0);
    fv.table.resize(kOutflowRows);
    fv.tableStep = window / (kOutflowRows - 1);
    fv.tableStart = fv.clock - window;
    double widest = 0.0;
    for (int k = 0; k < kOutflowRows; ++k) {
        const double t = ref.time + fv.tableStart + k * fv.tableStep;
        const double open = after.opening(t);
        widest = std::max(widest, open);
        fv.table[k] = glm::vec4(static_cast<float>(after.outflowRate(t)), static_cast<float>(open),
                                static_cast<float>(after.luminosity(t)), 0.0f);
    }
    fv.maxOpening = widest;

    // The widest the cone can be seen at each radius: gas there left a travel time ago,
    // seen a light-travel time (up to the camera's distance plus R) earlier still. The
    // billows reach 1.6 times the cone's half-width.
    fv.cone.assign(FlareView::kConeRows, 0.0f);
    const double front = after.front(now);
    for (int k = 0; k < FlareView::kConeRows; ++k) {
        const double R = R0 * std::pow(fv.extent / R0, static_cast<double>(k) / (FlareView::kConeRows - 1));
        if (R > front * 1.12 + 6.0) continue;
        const double latest = now - tde::windTravelTime(1.0, s.outflowSpeed, R0, R);
        const double earliest = latest - (s.lightDelay ? m_camera.r + R : 0.0) - 20.0;
        double open = 0.0;
        for (int j = 0; j <= 24; ++j) open = std::max(open, after.opening(earliest + (latest - earliest) * j / 24.0));
        fv.cone[k] = static_cast<float>(std::min(1.6 * open + 0.04, 1.5));
    }
    return fv;
}

void App::addEinsteinRing() {
    // Directly behind the hole as seen from the camera, at rest relative to the
    // local ZAMO: it falls almost straight in, staying lined up for a while.
    const double spin = m_settings.spin;
    const auto c = kerr::toCartesian(spin, m_camera.r, m_camera.theta, m_camera.phi);
    const double scale = -40.0 * sceneScale() / std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    const auto [r, theta, phi] = kerr::fromCartesian(spin, c[0] * scale, c[1] * scale, c[2] * scale);
    const double v[3] = {0.0, 0.0, 0.0};
    if (m_bodies.addFromZamo(spin, r, theta, phi, v, BodyKind::Star, "Ring star")) {
        m_bodies.bodies.back().temperature = 9000.0f;
        m_bodies.bodies.back().radius = 1.6f * static_cast<float>(sceneScale());
        m_bodies.logEvent("A star sits exactly behind the hole: its light forms an Einstein ring.");
    }
}

json::Value App::saveAppState() const {
    using json::Value;
    Value out = Value::object();
    out.set("flashTime", Value::number(m_flashTime));
    out.set("flashStrength", Value::number(m_flashStrength));
    out.set("flashFade", Value::number(m_flashFade));
    out.set("burstFuel", Value::number(m_burstFuel));
    out.set("lastAccreted", Value::number(m_lastAccreted));
    out.set("screenFlash", Value::number(m_screenFlash));
    Value flares = Value::array();
    for (const StarFlare& f : m_flares) {
        Value x = Value::object();
        x.set("direction", encodeVec3(f.direction));
        x.set("time", Value::number(f.time));
        x.set("strength", Value::number(f.strength));
        flares.push(std::move(x));
    }
    out.set("flares", std::move(flares));
    Value quakes = Value::array();
    for (const Quake& q : m_quakes) {
        Value x = Value::object();
        x.set("time", Value::number(q.time));
        x.set("amplitude", Value::number(q.amplitude));
        x.set("l", Value::number(static_cast<long long>(q.l)));
        x.set("direction", encodeVec3(q.direction));
        quakes.push(std::move(x));
    }
    out.set("quakes", std::move(quakes));
    std::ostringstream rng;
    rng << m_rng;
    out.set("random", Value::string(rng.str()));
    return out;
}

bool App::loadAppState(const json::Value& state, bool checkOnly, std::string& error) {
    // Each value keeps what it has if the file lacks it.
    double flashTime = m_flashTime, flashStrength = m_flashStrength, flashFade = m_flashFade;
    double burstFuel = 0.0, lastAccreted = 0.0;
    float screenFlash = 0.0f;
    std::vector<StarFlare> flares;
    std::vector<Quake> quakes;
    std::string random;
    auto read = [&](const char* key, auto& value) {
        const json::Value* x = state.find(key);
        return !x || x->get(value);
    };
    bool ok = state.isObject() && read("flashTime", flashTime) && read("flashStrength", flashStrength) &&
              read("flashFade", flashFade) && read("burstFuel", burstFuel) && read("lastAccreted", lastAccreted) &&
              read("screenFlash", screenFlash) && read("random", random);
    if (const json::Value* list = ok ? state.find("flares") : nullptr) {
        ok = list->isArray();
        for (size_t i = 0; ok && i < list->size(); ++i) {
            const json::Value& x = (*list)[i];
            StarFlare f{};
            const json::Value *time = x.find("time"), *strength = x.find("strength");
            ok = decodeVec3(x.find("direction"), f.direction) && time && time->get(f.time) && strength &&
                 strength->get(f.strength);
            flares.push_back(f);
        }
    }
    if (const json::Value* list = ok ? state.find("quakes") : nullptr) {
        ok = list->isArray();
        for (size_t i = 0; ok && i < list->size(); ++i) {
            const json::Value& x = (*list)[i];
            Quake q{};
            long long l = 2;
            const json::Value *time = x.find("time"), *amplitude = x.find("amplitude"), *degree = x.find("l");
            ok = time && time->get(q.time) && amplitude && amplitude->get(q.amplitude) && degree && degree->get(l) &&
                 l >= 2 && l <= 10 && decodeVec3(x.find("direction"), q.direction);
            q.l = static_cast<int>(l);
            quakes.push_back(q);
        }
    }
    std::mt19937 rng = m_rng;
    if (ok && !random.empty()) {
        std::istringstream in(random);
        in >> rng;
        ok = !in.fail();
    }
    if (!ok) {
        error = "the scene's flashes, flares or starquakes are damaged";
        return false;
    }
    if (checkOnly) return true;
    m_flashTime = flashTime;
    m_flashStrength = flashStrength;
    m_flashFade = flashFade;
    m_burstFuel = burstFuel;
    m_lastAccreted = lastAccreted;
    m_screenFlash = screenFlash;
    m_flares = std::move(flares);
    m_quakes = std::move(quakes);
    m_rng = rng;
    return true;
}

bool App::saveScene(const std::filesystem::path& path) {
    const std::string name = utf8Of(path.filename());
    if (m_merger.active()) {
        m_ui.sceneStatus = "The neutron star merger runs a simulation of its own and cannot be saved; leave it first.";
        return false;
    }
    // Debris added since the last step goes to the gas now, so the file holds all of it.
    for (const auto& batch : m_bodies.pendingDebris) m_particles.spawn(batch.particles, m_bodies.time, batch.existedBefore);
    m_bodies.pendingDebris.clear();

    scene::Blob blob;
    json::Value root = scene::save(m_settings, m_camera, m_bodies, blob);
    root.set("gas", m_particles.save(blob));
    root.set("app", saveAppState());
    root.set("keyframes", keyframes::toJson(m_track.keys));
    std::string error;
    if (!scene::writeFiles(path, std::move(root), blob, error)) {
        m_ui.sceneStatus = "Could not save the scene: " + error;
        std::fprintf(stderr, "%s\n", m_ui.sceneStatus.c_str());
        return false;
    }
    m_scenePath = path;
    m_ui.sceneName = name;
    m_ui.sceneStatus = "Saved " + name + " (with " + name + ".bin beside it: the bodies' past and the gas).";
    std::printf("Saved scene %s\n", utf8Of(path).c_str());
    return true;
}

bool App::loadScene(const std::filesystem::path& path) {
    const std::string name = utf8Of(path.filename());
    json::Value root;
    scene::Blob blob;
    std::string error;
    auto fail = [&]() {
        m_ui.sceneStatus = "Could not open " + name + ": " + error;
        std::fprintf(stderr, "%s\n", m_ui.sceneStatus.c_str());
        return false;
    };
    if (!scene::readFiles(path, root, blob, error)) return fail();
    // Everything is checked before anything changes; the bodies, settings and camera
    // are only replaced if they all read correctly.
    const json::Value* gas = root.find("gas");
    const json::Value* app = root.find("app");
    if (gas && !m_particles.load(*gas, blob, error, true)) return fail();
    if (app && !loadAppState(*app, true, error)) return fail();
    if (!scene::load(root, blob, m_settings, m_camera, m_bodies, error)) return fail();

    if (m_merger.active()) m_merger.stop();
    clearTimeline();
    if (gas) m_particles.load(*gas, blob, error);
    else m_particles.clear();
    if (app) loadAppState(*app, false, error);
    std::vector<keyframes::Keyframe> keys;
    if (const json::Value* list = root.find("keyframes"); list && keyframes::fromJson(*list, keys)) {
        m_track.keys = std::move(keys);
        m_track.playing = false;
        m_track.rebuild();
    }
    // The loaded spin and scale are the ones everything is already in.
    m_lastSpin = m_settings.spin;
    m_lastScale = sceneScale();
    // Predictions are for the bodies that were here before.
    if (m_jointJob.valid()) m_jointJob.wait();
    m_jointJob = {};
    m_jointSlots.clear();
    m_jointPlannerKey.clear();
    m_predictions.assign(History::kSlots, trajectory::Path{});
    m_predictNext = 0;
    m_streamQ = m_streamRim = -1.0;
    m_ui.compactnessOrigin = -1;
    m_ui.inspectValid = false;
    m_scenePath = path;
    m_ui.sceneName = name;
    char text[160];
    std::snprintf(text, sizeof(text), "Opened %s: %zu bodies, t = %.6g M.", name.c_str(), m_bodies.bodies.size(),
                  m_bodies.time);
    m_ui.sceneStatus = text;
    m_ui.cameraStatus.clear();
    m_ui.centralStatus.clear();
    std::printf("%s\n", text);
    return true;
}

void App::saveSceneInteractive(bool choose) {
    if (!choose && !m_scenePath.empty()) {
        saveScene(m_scenePath);
        return;
    }
    const std::filesystem::path folder = m_scenePath.empty() ? std::filesystem::path("scenes") : m_scenePath.parent_path();
    const std::filesystem::path name = m_scenePath.empty()
        ? std::filesystem::path(timestampedPath("scenes", scene::kExtension)).filename()
        : m_scenePath.filename();
    if (auto path = dialog::saveScene(m_window, folder, name)) saveScene(*path);
}

void App::openSceneInteractive() {
    const std::filesystem::path folder = m_scenePath.empty() ? std::filesystem::path("scenes") : m_scenePath.parent_path();
    if (auto path = dialog::openScene(m_window, folder)) loadScene(*path);
}

void App::simulateFrame(double simDt, float realDt) {
    syncCentral(simDt, realDt);
    m_bodies.setFrameTime(simDt, m_settings.spin);
    advanceSimulation(simDt);
    syncRide(simDt);
    if (m_settings.flareCamera && m_camera.motion == CameraMotion::Hover) {
        const double target = flareCameraDistance();
        m_camera.r += (target - m_camera.r) * (1.0 - std::exp(-realDt / 1.5));
    }
    if (m_camera.lookAhead && m_camera.motion != CameraMotion::Hover) m_camera.faceVelocity(m_settings.spin);
}

std::uint64_t App::simulationFingerprint() const {
    std::uint64_t hash = 14695981039346656037ull;  // FNV-1a
    auto mix = [&](const auto& value) {
        const unsigned char* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (size_t i = 0; i < sizeof(value); ++i) {
            hash ^= bytes[i];
            hash *= 1099511628211ull;
        }
    };
    mix(m_bodies.time);
    mix(m_bodies.bodies.size());
    for (const Body& b : m_bodies.bodies) {
        mix(b.slot);
        mix(b.kind);
        mix(b.mass);
        mix(b.radius);
        mix(b.temperature);
        mix(b.brightness);
        mix(b.spinRate);
        mix(b.tidalRadius);
        mix(b.constants);
        mix(b.state);
        mix(b.shape);
    }
    mix(m_bodies.departed.size());
    mix(m_bodies.tidalEvents.size());
    mix(m_bodies.pendingDebris.size());
    mix(m_bodies.impacts.size());
    mix(m_particles.edits());
    mix(m_flashTime);
    mix(m_burstFuel);
    mix(m_quakes.size());
    mix(m_flares.size());
    mix(m_camera.motion);
    if (m_camera.motion != CameraMotion::Hover) {
        mix(m_camera.path);
        mix(m_camera.geodesic);
        mix(m_camera.rideSlot);
    }
    return hash;
}

void App::addCheckpoint(bool pinned) {
    // Debris waiting for the gas goes to it now, so the checkpoint holds all of it.
    for (const auto& batch : m_bodies.pendingDebris) m_particles.spawn(batch.particles, m_bodies.time, batch.existedBefore);
    m_bodies.pendingDebris.clear();
    Timeline::Checkpoint c;
    c.pinned = pinned;
    json::Value root = scene::save(m_settings, m_camera, m_bodies, c.blob);
    root.set("app", saveAppState());
    c.scene = json::dump(root);
    ParticleSystem::GpuCopy copy = m_particles.copyOnGpu();
    if (copy.allocated) {
        c.gpuBytes = copy.bytes;
        size_t slot = 0;
        while (slot < m_gasCopies.size() && m_gasCopies[slot].allocated) ++slot;
        if (slot == m_gasCopies.size()) m_gasCopies.emplace_back();
        m_gasCopies[slot] = copy;
        c.gpu = static_cast<int>(slot);
    }
    if (m_timeline.empty()) m_timeline.setSettings(scene::simulationSettings(m_settings));
    const int replaced = m_timeline.addCheckpoint(std::move(c), m_bodies.time);
    if (replaced >= 0) releaseGasCopies({replaced});
}

void App::releaseGasCopies(const std::vector<int>& slots) {
    for (int slot : slots) {
        if (slot < 0 || slot >= static_cast<int>(m_gasCopies.size())) continue;
        ParticleSystem::GpuCopy& copy = m_gasCopies[static_cast<size_t>(slot)];
        if (copy.buffer) glDeleteBuffers(1, &copy.buffer);
        copy = ParticleSystem::GpuCopy{};
    }
}

void App::clearTimeline() {
    for (ParticleSystem::GpuCopy& copy : m_gasCopies) {
        if (copy.buffer) glDeleteBuffers(1, &copy.buffer);
    }
    m_gasCopies.clear();
    m_timeline.clear();
    m_seekTarget = -1;
    m_reverse = false;
}

bool App::restoreCheckpoint(const Timeline::Checkpoint& c) {
    json::Value root;
    std::string error;
    Settings settings = m_settings;
    Camera camera = m_camera;
    if (!json::parse(c.scene, root, error) || !scene::load(root, c.blob, settings, camera, m_bodies, error)) {
        m_ui.cameraStatus = "The timeline could not go back: " + error;
        return false;
    }
    // What moves comes back as it was; how it looks stays as you have it now. A camera on
    // a path (falling, orbiting, riding) is part of what moves.
    scene::applySimulationSettings(scene::simulationSettings(settings), m_settings);
    m_settings.spin = settings.spin;
    if (camera.motion != CameraMotion::Hover || m_camera.motion != CameraMotion::Hover) m_camera = camera;
    if (c.gpu >= 0) m_particles.restoreFromGpu(m_gasCopies[static_cast<size_t>(c.gpu)]);
    else m_particles.restoreFromGpu(ParticleSystem::GpuCopy{});
    if (const json::Value* app = root.find("app")) loadAppState(*app, false, error);
    m_lastSpin = m_settings.spin;
    m_lastScale = sceneScale();
    if (m_jointJob.valid()) m_jointJob.wait();
    m_jointJob = {};
    m_jointSlots.clear();
    m_jointPlannerKey.clear();
    m_predictions.assign(History::kSlots, trajectory::Path{});
    m_predictNext = 0;
    m_timeline.setPosition(c.step);
    return true;
}

bool App::seekTowards(size_t target, double seconds) {
    target = std::min(target, m_timeline.length());
    size_t position = m_timeline.position();
    const Timeline::Checkpoint& c = m_timeline.checkpointFor(target);
    // Going back, or a checkpoint nearer than where the scene is: start from it.
    if (target < position || c.step > position) {
        if (!restoreCheckpoint(c)) return true;
        position = c.step;
    }
    const auto start = std::chrono::steady_clock::now();
    while (position < target) {
        const Timeline::Step& step = m_timeline.step(position);
        if (!step.settings.empty()) scene::applySimulationSettings(step.settings, m_settings);
        simulateFrame(step.dt, step.realDt);
        m_timeline.setPosition(++position);
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() > seconds) break;
    }
    m_timeline.setSettings(scene::simulationSettings(m_settings));
    return position == target;
}

void App::timelineFrame(float realDt, double liveDt) {
    // What did the user change since the last frame?
    const std::string settings = scene::simulationSettings(m_settings);
    const bool started = !m_timeline.empty();
    const bool acted = started && simulationFingerprint() != m_fingerprint;
    const bool retuned = started && settings != m_timeline.settings();
    if ((acted || retuned) && !m_timeline.atEnd()) {
        // Changed while looking at the past: a new future starts here.
        releaseGasCopies(m_timeline.truncate());
        m_seekTarget = -1;
        m_reverse = false;
        m_ui.cameraStatus = "Timeline: a new future starts here; what came after this moment was dropped.";
    }
    if (!started || acted) addCheckpoint(true);

    bool stepped = false;
    if (m_seekTarget >= 0) {
        // Scrubbing: replay toward the moment asked for, a little each frame.
        stepped = true;
        if (seekTowards(static_cast<size_t>(m_seekTarget), kSeekBudget)) m_seekTarget = -1;
    } else if (!m_timeline.atEnd() || (m_reverse && !m_settings.paused)) {
        // In the past: play on, or back, one recorded frame per frame.
        if (!m_settings.paused) {
            const size_t position = m_timeline.position();
            if (m_reverse) {
                if (position == 0) {
                    m_settings.paused = true;
                    m_reverse = false;
                } else {
                    m_seekTarget = static_cast<long long>(position - 1);
                    if (seekTowards(position - 1, kSeekBudget)) m_seekTarget = -1;
                    stepped = true;
                }
            } else {
                stepped = true;
                seekTowards(position + 1, kSeekBudget);
            }
        }
    } else {
        // Live: step the scene and record the step.
        m_reverse = false;
        simulateFrame(liveDt, realDt);
        stepped = true;
        if (liveDt > 0.0) {
            Timeline::Step step;
            step.dt = liveDt;
            step.realDt = realDt;
            if (retuned) step.settings = settings;
            m_timeline.addStep(std::move(step), m_bodies.time);
            m_timeline.setSettings(scene::simulationSettings(m_settings));
            const double tScale = sceneScale() * std::sqrt(sceneScale());
            if (m_bodies.time - m_timeline.lastCheckpointTime() >= kCheckpointInterval * tScale) addCheckpoint(false);
            releaseGasCopies(m_timeline.thin(kTimelineMemory, kTimelineGpuMemory));
        }
    }
    // Paused, or holding still in the past: what each frame does anyway (the gas grids are
    // rebuilt after a restore, flashes fade).
    if (!stepped) simulateFrame(0.0, realDt);
    m_fingerprint = simulationFingerprint();
}

UiContext::TimelineView App::timelineView() const {
    UiContext::TimelineView v;
    v.available = m_timelineOn && !m_timeline.empty() && !m_merger.active();
    v.start = m_timeline.start();
    v.end = m_timeline.end();
    v.time = m_timeline.time();
    v.live = m_timeline.atEnd() && m_seekTarget < 0 && !m_reverse;
    v.playing = !m_settings.paused;
    v.reverse = m_reverse;
    v.seeking = m_seekTarget >= 0 && !m_reverse;
    v.checkpoints = m_timeline.checkpointCount();
    v.megabytes = static_cast<double>(m_timeline.bytes()) / (1 << 20);
    v.gpuMegabytes = static_cast<double>(m_timeline.gpuBytes()) / (1 << 20);
    return v;
}

void App::timelineCommand(int command) {
    if (!m_timelineOn || m_timeline.empty()) return;
    const size_t position = m_timeline.position();
    switch (command) {
        case 1:  // play or pause (forward)
            if (m_reverse) {
                m_reverse = false;
                m_settings.paused = false;
            } else {
                m_settings.paused = !m_settings.paused;
            }
            break;
        case 2:  // play backwards (or stop doing so)
            m_reverse = !m_reverse || m_settings.paused;
            m_settings.paused = !m_reverse;
            m_seekTarget = -1;
            break;
        case 3:  // back to the present
            m_reverse = false;
            m_seekTarget = static_cast<long long>(m_timeline.length());
            break;
        case 4:  // to the start
            m_reverse = false;
            m_settings.paused = true;
            m_seekTarget = 0;
            break;
        case 5:  // one frame back
            m_reverse = false;
            m_settings.paused = true;
            m_seekTarget = static_cast<long long>(position > 0 ? position - 1 : 0);
            break;
        case 6:  // one frame on
            m_reverse = false;
            m_settings.paused = true;
            m_seekTarget = static_cast<long long>(std::min(position + 1, m_timeline.length()));
            break;
        default:
            break;
    }
}

std::vector<SceneSound::Voice> App::soundVoices(double timeRate, double& rumble, double& waveHz) const {
    std::vector<SceneSound::Voice> voices;
    rumble = 0.0;
    waveHz = 0.0;
    const double a = m_settings.spin, scale = sceneScale();
    if (m_camera.insideHorizon(a)) return voices;
    // The camera: where it is, which way is right, and its clock's rate (a hovering
    // camera's ZAMO lapse).
    const auto pc = kerr::toCartesian(a, m_camera.r, m_camera.theta, m_camera.phi);
    const glm::dvec3 camera(pc[0], pc[1], pc[2]);
    const glm::dmat3 view = glm::mat3_cast(keyframes::fromCamera(m_camera, a).orientation);
    const glm::dvec3 right = view[0];
    const kerr::Metric mc = kerr::metricAt(a, m_camera.r, m_camera.theta);
    const double lapseCamera = std::sqrt(std::max(mc.sigma * mc.delta / mc.A, 1e-12));
    double loudestWave = 0.0;
    for (const Body& b : m_bodies.bodies) {
        if (b.inside) continue;
        const auto pb = kerr::toCartesian(a, b.state.r, b.state.theta, b.state.phi);
        const glm::dvec3 body(pb[0], pb[1], pb[2]);
        const glm::dvec3 toCamera = camera - body;
        const double distance = std::max(glm::length(toCamera), 1e-6);
        const glm::dvec3 n = toCamera / distance;
        // g: the body's clock against the camera's (lapses), times the Doppler shift of its
        // velocity relative to its ZAMO along the line of sight (light bending left out).
        double v[3] = {0.0, 0.0, 0.0};
        if (b.state.r > kerr::horizonRadius(a) + 0.05) kerr::zamoVelocity(a, b.constants, b.state, v);
        const kerr::LocalAxes axes = kerr::localAxes(a, b.state.r, b.state.theta, b.state.phi);
        glm::dvec3 vel(0.0);
        for (int i = 0; i < 3; ++i) vel[i] = v[0] * axes.r[i] + v[1] * axes.theta[i] + v[2] * axes.phi[i];
        const double speed2 = std::min(glm::dot(vel, vel), 0.999999);
        const double doppler = std::sqrt(1.0 - speed2) / (1.0 - glm::dot(vel, n));
        const kerr::Metric mb = kerr::metricAt(a, b.state.r, b.state.theta);
        const double lapseBody = std::sqrt(std::max(mb.sigma * mb.delta / mb.A, 0.0));
        const double g = std::clamp(lapseBody / lapseCamera * doppler, 0.05, 20.0);
        const double near = 1.0 / (1.0 + (distance / (25.0 * scale)) * (distance / (25.0 * scale)));
        const double pan = std::clamp(glm::dot(-n, right), -1.0, 1.0);
        const int id = 16 * b.slot;
        if (b.kind == BodyKind::Pulsar && m_settings.soundPulsars) {
            // Two beams: two ticks a turn, at the shifted rate, louder by g^3.
            const double rate = std::abs(b.spinRate) * timeRate / kerr::kPi * g;
            voices.push_back({id, SceneSound::Clicks, std::min(2600.0 * g, 9000.0),
                              std::min(0.5 * near * g * g * g, 1.0), pan, std::min(rate, 400.0)});
        } else if (m_settings.soundBodies && b.kind != BodyKind::MiniBlackHole) {
            double pitch = 220.0;  // Hz, before the shift
            switch (b.kind) {
                case BodyKind::Star:
                case BodyKind::TidalStar: pitch = 220.0 * std::pow(std::max(b.temperature, 1000.0f) / 5800.0, 0.35); break;
                case BodyKind::RockyPlanet: pitch = 131.0; break;
                case BodyKind::GasGiant: pitch = 98.0; break;
                case BodyKind::Moon: pitch = 392.0; break;
                case BodyKind::WhiteDwarf: pitch = 523.0; break;
                case BodyKind::NeutronStar: pitch = 659.0; break;
                default: break;
            }
            voices.push_back({id, SceneSound::Tone, std::clamp(pitch * g, 25.0, 8000.0),
                              std::min(0.25 * near * g * g, 0.7), pan, 0.0});
        }
        // Gravitational waves at twice the orbital frequency, from bodies heavy enough to
        // send out any worth hearing; raised in pitch into hearing.
        if (m_settings.soundWaves && b.mass >= 1e-3 && m_settings.central.hasHorizon()) {
            const kerr::Derivatives d = kerr::derivatives(a, b.constants, b.state);
            const double omega = std::abs(d.dphi / d.dt);  // radians per M
            const double heard = omega / kerr::kPi * timeRate;  // cycles per real second
            const double amplitude = std::min(1.0, 30.0 * b.mass) * std::pow(std::min(omega * 14.7, 1.0), 2.0 / 3.0);
            voices.push_back({id + 1, SceneSound::Wave, std::clamp(heard * 400.0, 30.0, 1500.0),
                              0.25 * amplitude * near, pan, 0.0});
            if (amplitude > loudestWave) {
                loudestWave = amplitude;
                waveHz = omega / kerr::kPi / central::timeSeconds(1.0, m_settings.central.massSolar);
            }
        }
    }
    // The tide at the camera grows as 1/r^3: a roar, heard from about the ISCO in.
    if (m_settings.soundRumble && m_settings.central.hasHorizon()) {
        const double x = 6.0 / std::max(m_camera.r, 1e-3);
        rumble = std::clamp(0.35 * x * x * x, 0.0, 1.0);
    }
    return voices;
}

EhtView::Options App::ehtOptions() const {
    EhtView::Options o;
    const UiState& ui = m_ui;
    o.target = std::clamp(ui.ehtTarget, 0, 2);
    if (o.target == EhtView::ThisHole) {
        o.spin = m_settings.central.hasHorizon() ? m_settings.spin : 0.0;
        o.inclinationDeg = std::clamp(m_camera.theta * 180.0 / kerr::kPi, 1.0, 179.0);
        o.massSolar = m_settings.central.massSolar;
    } else {
        o.spin = ui.ehtSpin;
        o.inclinationDeg = ui.ehtInclination[o.target];
    }
    const double frequencies[3] = {86.0, 230.0, 345.0};
    o.frequencyGHz = frequencies[std::clamp(ui.ehtFrequency, 0, 2)];
    o.temperature = ui.ehtTemperature[o.target];
    o.beta = ui.ehtBeta;
    o.height = ui.ehtHeight[o.target];
    const int sizes[3] = {96, 160, 256};
    o.size = m_ehtSizeOverride > 0 ? m_ehtSizeOverride : sizes[std::clamp(ui.ehtResolution, 0, 2)];
    o.beamUas = ui.ehtBeam;
    o.curve = ui.ehtCurve;
    return o;
}

void App::addKeyframe() {
    if (m_merger.active()) return;
    if (m_camera.insideHorizon(m_settings.spin)) {
        m_track.status = "Keyframes are places a camera can hover: not inside the horizon.";
        return;
    }
    keyframes::Keyframe k = keyframes::fromCamera(m_camera, m_settings.spin);
    // A first guess at the time to get here: about 3 s, longer for a long way.
    if (!m_track.keys.empty()) {
        const double far = glm::length(k.position - m_track.keys.back().position) / (10.0 * sceneScale());
        m_track.keys.back().duration = std::clamp(3.0 * std::sqrt(std::max(far, 0.1)), 2.0, 12.0);
    }
    m_track.keys.push_back(k);
    m_track.rebuild();
    char text[96];
    std::snprintf(text, sizeof(text), "Keyframe %zu added.", m_track.keys.size());
    m_track.status = text;
}

void App::playKeyframes(float dt, double timeRate) {
    const bool filming = m_recording.active && m_recording.cameraPath == RecordKeyframes;
    if ((!m_track.playing && !filming) || m_track.path.empty() || m_merger.active()) {
        m_track.playing = m_track.playing && !m_track.path.empty();
        m_camera.drift = {};
        return;
    }
    if (filming) {
        m_track.time = static_cast<double>(m_recording.frame) / m_recording.fps;
    } else {
        m_track.time += dt;
        if (m_track.time > m_track.path.duration()) {
            m_track.playing = false;
            m_camera.drift = {};
            m_track.status = "The path is over.";
            return;
        }
    }
    const keyframes::Pose pose = m_track.path.at(std::min(m_track.time, m_track.path.duration()));
    if (m_camera.motion != CameraMotion::Hover) setCameraMotion(CameraMotion::Hover);
    m_settings.flareCamera = false;  // the path sets the distance
    keyframes::applyPose(pose, m_settings.spin, m_camera.minDistance(m_settings.spin), m_camera);
    if (m_track.aim) m_camera.yaw = m_camera.pitch = 0.0f;  // facing the center, as yaw = pitch = 0 does
    m_camera.drift = {};
    if (m_track.physical) {
        if (timeRate <= 0.0) {
            m_track.status = "Time is paused: nothing can move for real. The camera hovers at each point.";
        } else {
            // The path's speed in the simulation's time (M per M: in c).
            bool capped = false;
            const glm::dvec3 v = keyframes::zamoVelocity(pose.position, pose.velocity / timeRate, m_settings.spin,
                                                         0.99, capped);
            m_camera.drift = {v.x, v.y, v.z};
            char text[160];
            std::snprintf(text, sizeof(text), "Moving at %.3f c%s.", glm::length(v),
                          capped ? " (capped: the path is faster than light at this pace of time)" : "");
            m_track.status = text;
        }
    }
}

void App::advanceSimulation(double dt) {
    const double spin = m_settings.spin;
    m_bodies.update(spin, dt, m_settings.bodyGravity);
    // The Penrose process: once a negative-energy fragment has fallen in, the hole has lost
    // mass m1 E1 and angular momentum m1 L1, and turns more slowly (src/penrose.h).
    for (penrose::Split& split : m_bodies.penroseSplits) {
        if (!split.absorbed || split.applied) continue;
        split.applied = true;
        if (!m_settings.penroseSpinDown || !m_settings.central.hasHorizon()) continue;
        const double after = std::clamp(penrose::spinAfter(spin, split), 0.0, 0.999);
        m_settings.central.chi = after;
        char text[300];
        std::snprintf(text, sizeof(text),
                      "The negative-energy fragment fell in: the hole lost %.2g of its mass and %.2g of its angular "
                      "momentum. Its spin fell from %.4f to %.4f.",
                      -split.m1 * split.E1, -split.m1 * split.L1 / std::max(spin, 1e-9), spin, after);
        m_bodies.logEvent(text);
    }
    for (const auto& batch : m_bodies.pendingDebris) m_particles.spawn(batch.particles, m_bodies.time, batch.existedBefore);
    m_bodies.pendingDebris.clear();
    {
        // The unbound half of a star torn apart on a parabolic orbit flies far out: give the
        // gas grid room for it.
        bool far = false;
        for (const tde::Event& e : m_bodies.tidalEvents) far = far || e.unboundFraction > 0.15;
        m_particles.setGridOuterRadius(far && m_settings.central.hasHorizon() ? kFarDebrisRadius
                                                                              : ParticleSystem::kGridRMax);
    }
    GasParams gas;
    const CentralView view = centralView();
    gas.scale = sceneScale();
    gas.surface = view.surface;
    gas.surfaceShift = view.azimuthShift;
    gas.surfacePhase = view.phase;
    gas.impactFadeTime = view.impactFade;
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
    if (m_camera.landed) {
        // Touchdown on a star or planet.
        m_camera.landed = false;
        const central::Object& o = m_settings.central;
        const double v = m_camera.landingSpeed;
        const double g = central::accelerationSI(central::surfaceGravity(o.radius), o.massSolar) / 9.80665;
        char text[400];
        std::snprintf(text, sizeof(text),
                      "You landed on the surface at %.3g c (gamma %.3g). Surface gravity: %.3g g. Clocks here run at "
                      "%.4g%% of the rate far away.",
                      v, 1.0 / std::sqrt(std::max(1.0 - v * v, 1e-12)), g, 100.0 * central::surfaceClockRate(o.radius));
        m_ui.cameraStatus = text;
    }
    if (m_camera.insideHorizon(spin) && !wasInside && m_camera.motion != CameraMotion::Hover) {
        m_ui.cameraStatus = "You crossed the event horizon. Nothing marked the moment, but there is no way "
                            "back: from here every path, even light's, leads inward.";
    }
    if (m_camera.ended && !wasEnded) {
        m_ui.cameraStatus = spin > 1.0
            ? "The fall ends beside the ring singularity, where the curvature grows without limit. Press R to start "
              "over."
            : spin > 0.01
            ? "The fall ends just above the inner horizon. There the Kerr interior stops being trustworthy "
              "(infalling light piles up without limit), so the simulation stops. Press R to start over."
            : "The fall ends just above the singularity, where tides crush everything. Press R to start over.";
        m_settings.paused = true;
    }
}

bool App::startRecording(const std::string& path, int width, int height, int fps, float quality, int cameraPath,
                         double seconds) {
    if (cameraPath == RecordKeyframes) {
        if (m_track.keys.size() < 2) {
            m_ui.videoStatus = "Add at least two camera keyframes (K) to film along them.";
            std::fprintf(stderr, "%s\n", m_ui.videoStatus.c_str());
            return false;
        }
        seconds = m_track.path.duration();
        m_track.playing = false;
    }
    ensureParent(path);
    std::string error;
    const int audioRate = m_merger.active() ? (m_merger.sound ? MergerScene::kSampleRate : 0)
                                            : (m_settings.sound ? SceneSound::kSampleRate : 0);
    if (!m_recording.encoder.open(path, width, height, fps, error, audioRate)) {
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
    m_recording.soundSquares = 0.0;
    m_recording.soundSamples = 0;
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
    if (m_merger.active() || m_settings.sound) {
        const std::vector<int16_t>& sound = m_merger.active() ? m_merger.lastAudio() : m_sound.lastMono();
        for (int16_t x : sound) m_recording.soundSquares += static_cast<double>(x) * x;
        m_recording.soundSamples += static_cast<long long>(sound.size());
        if (!m_recording.encoder.writeAudio(sound.data(), sound.size(), error)) std::fprintf(stderr, "%s\n", error.c_str());
    }
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
        int n = std::snprintf(text, sizeof(text), "Saved %s (%d frames, %.1f s)", m_recording.path.c_str(),
                              m_recording.frame, static_cast<double>(m_recording.frame) / m_recording.fps);
        if (m_recording.soundSamples > 0 && n > 0) {
            // How loud the sound track is (rms, as a fraction of full scale).
            std::snprintf(text + n, sizeof(text) - n, ", with sound (%.1f s, rms %.3f)",
                          static_cast<double>(m_recording.soundSamples) / SceneSound::kSampleRate,
                          std::sqrt(m_recording.soundSquares / m_recording.soundSamples) / 32768.0);
        }
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

        // A scene dropped onto the window (or its .bin sidecar) is opened.
        if (!m_inbox.dropped.empty()) {
            std::filesystem::path path = pathFromUtf8(m_inbox.dropped.front());
            m_inbox.dropped.clear();
            if (path.extension() == ".bin") path.replace_extension();
            if (path.extension() == scene::kExtension) {
                loadScene(path);
            } else {
                m_ui.sceneStatus = "Only scene files (.bhscene) can be dropped here.";
            }
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        if (m_merger.active()) m_camera.scale = 1e15;  // zoom out as far as the kilonova
        if (!headless) handleInput(dt);

        // Simulation. Time runs in scene units: the dynamical time grows as scale^1.5.
        // The neutron star merger runs a scene of its own (src/merger_scene.h).
        const bool merging = m_merger.active();
        if (merging) {
            m_camera.scale = 1e15;
            const int budget = headless || m_recording.active ? 600 : 14;
            m_merger.update(dt, m_settings.paused || screenshotMode, m_camera, budget, !headless && !m_recording.active);
        }
        const bool advance = !merging && !m_settings.paused && !screenshotMode;
        const double tScale = sceneScale() * std::sqrt(sceneScale());
        const double simDt = advance ? dt * m_settings.timeScale * tScale * approachSlowdown() * flarePace() : 0.0;
        if (!merging) playKeyframes(dt, dt > 0.0f ? simDt / dt : 0.0);
        if (m_timelineOn && !merging && !screenshotMode) timelineFrame(dt, simDt);
        else simulateFrame(simDt, dt);
        {
            // The scene's sound, carried by light (the merger makes its own).
            const bool on = m_settings.sound && !merging;
            double rumble = 0.0;
            std::vector<SceneSound::Voice> voices;
            if (on) voices = soundVoices(dt > 0.0f ? simDt / dt : 0.0, rumble, m_waveHz);
            m_sound.update(voices, rumble, dt, !headless && !m_recording.active, on, m_settings.soundVolume);
        }
        updatePredictions();

        if (merging) {
            // The merger's camera follows the action by itself.
        } else if (m_recording.active && m_camera.motion == CameraMotion::Hover) {
            const double progress = static_cast<double>(m_recording.frame) / m_recording.total;
            if (m_recording.cameraPath == RecordOrbit) {
                m_camera.phi = m_recording.startPhi + 2.0 * kerr::kPi * progress;
            } else if (m_recording.cameraPath == RecordPushIn) {
                m_camera.r = std::max(m_recording.startR * std::pow(0.3, progress),
                                      m_camera.minDistance(m_settings.spin) + 0.7 * sceneScale());
            }
        } else if (m_settings.autoOrbit && !headless && m_camera.motion == CameraMotion::Hover && !m_track.playing) {
            m_camera.orbit(m_settings.autoOrbitSpeed * dt, 0.0);
        }
        if (m_settings.vsync != m_vsync && !headless && !m_recording.active) {
            m_vsync = m_settings.vsync;
            glfwSwapInterval(m_vsync ? 1 : 0);
        }
        if (now - lastReloadCheck > 0.5) {
            m_renderer.hotReload();
            m_particles.hotReload();
            m_merger.hotReload();
            lastReloadCheck = now;
        }

        // Dynamic resolution: the traced resolution follows a frame-time target,
        // scale *= (target / measured)^(1/2) (the cost goes as the pixel count), smoothed and
        // stepped by 1/40 so the targets are not reallocated every frame.
        if (m_settings.dynamicResolution && !headless && !m_recording.active) {
            if (m_dynamicScale <= 0.0f) m_dynamicScale = m_settings.renderScale;
            const float measured = m_renderer.gpuMilliseconds();
            if (measured > 0.5f) {
                const float wanted =
                    std::clamp(m_settings.renderScale * std::sqrt(m_settings.targetMs / measured), 0.4f, 1.0f);
                m_dynamicScale += (wanted - m_dynamicScale) * 0.1f;
                const float stepped = std::round(m_dynamicScale * 40.0f) / 40.0f;
                if (std::abs(stepped - m_settings.renderScale) >= 0.025f) m_settings.renderScale = stepped;
            }
        } else {
            m_dynamicScale = 0.0f;
        }

        // Render
        const int outW = m_recording.active ? m_recording.width : fbw;
        const int outH = m_recording.active ? m_recording.height : fbh;
        if (merging) {
            m_merger.render(m_renderer, m_settings, m_camera, outW, outH,
                            m_recording.active ? m_recording.quality : m_settings.renderScale);
        } else {
            m_renderer.uploadHistory(m_bodies.history);
            buildPaths();
            FrameInputs in;
            in.width = m_recording.active ? m_recording.width : fbw;
            in.height = m_recording.active ? m_recording.height : fbh;
            if (m_recording.active) in.renderScale = m_recording.quality;
            GpuBodies gpuBodies;
            {
                // A traced pixel's angle, for drawing bodies too small to resolve as points.
                const double traceHeight =
                    std::max(1.0, in.height * static_cast<double>(in.renderScale > 0.0f ? in.renderScale
                                                                                        : m_settings.renderScale));
                const double pixelAngle = 2.0 * std::tan(0.5 * m_camera.fovDeg * kerr::kPi / 180.0) / traceHeight;
                const auto cam = kerr::toCartesian(m_settings.spin, m_camera.r, m_camera.theta, m_camera.phi);
                m_bodies.fillGpu(m_settings.spin, m_settings.lightDelay, m_camera.r, glm::dvec3(cam[0], cam[1], cam[2]),
                                 pixelAngle, gpuBodies);
            }
            const CentralView view = centralView();
            const FlareView flare = flareView();
            in.flare = &flare;
            in.settings = &m_settings;
            in.camera = &m_camera;
            in.bodies = &gpuBodies;
            in.particles = &m_particles;
            in.central = &view;
            in.skipBody = ridingIndex();
            in.spin = m_settings.spin;
            in.scale = sceneScale();
            in.reach = m_bodies.reach();
            in.exposureBoost = 1.0f + 6.0f * m_screenFlash;
            in.time = m_bodies.time;
            m_renderer.render(in);
        }
        if (m_recording.active) recordFrame();
        if (!headless) m_renderer.present(fbw, fbh);

        UiContext ui{m_settings, m_camera, m_bodies, m_renderer, m_particles, m_ui, m_fps, m_gpuName.c_str()};
        ui.recording = m_recording.active;
        ui.recordFrame = m_recording.frame;
        ui.recordTotal = m_recording.total;
        ui.predictions = &m_predictions;
        ui.slowdown = approachSlowdown();
        ui.pace = flarePace();
        ui.scale = sceneScale();
        ui.merger = &m_merger;
        ui.timeline = timelineView();
        ui.track = &m_track;
        ui.waveHz = m_waveHz;
        if (m_ui.ehtOpen && !merging) {
            m_eht.update(ehtOptions());
            if (screenshotMode) m_eht.finish();  // the picture must hold the finished image
            const bool panel = !headless || options.screenshotUi;
            const float em = ImGui::GetFontSize();
            const bool bar = panel && m_ui.showPanel && m_ui.showTimeline && ui.timeline.available;
            m_eht.draw(panel && m_ui.showPanel ? m_ui.panelRight + 0.5f * em : 0.0f, bar ? 5.6f * em : 0.0f);
        } else if (m_eht.busy()) {
            m_eht.shutdown();
        }
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
        if (ui.sceneFileRequest == 1 || ui.sceneFileRequest == 2) saveSceneInteractive(ui.sceneFileRequest == 2);
        if (ui.sceneFileRequest == 3) openSceneInteractive();
        if (ui.cameraMotionRequested >= 0) setCameraMotion(static_cast<CameraMotion>(ui.cameraMotionRequested));
        if (ui.rideRequested >= 0) rideWith(ui.rideRequested);
        if (ui.sceneRequested >= 0) applyScene(ui.sceneRequested);
        if (ui.einsteinRingRequested) addEinsteinRing();
        if (ui.centralRequested >= 0) applyCentral(ui.centralRequested);
        if (ui.quakeRequested) starquake();
        if (ui.timelineSeekRequested && m_timelineOn && !m_timeline.empty()) {
            m_reverse = false;
            m_settings.paused = true;
            m_seekTarget = static_cast<long long>(m_timeline.positionAt(ui.timelineSeekTime));
        }
        if (ui.timelineCommand) timelineCommand(ui.timelineCommand);
        if (ui.keyframeRequest == 1) addKeyframe();
        if (ui.keyframeRequest == 2) {
            m_track.playing = !m_track.playing && m_track.keys.size() >= 2;
            m_track.time = 0.0;
            m_track.status.clear();
            if (!m_track.playing) m_camera.drift = {};
        }
        if (ui.keyframeGo >= 0 && ui.keyframeGo < static_cast<int>(m_track.keys.size())) {
            m_track.playing = false;
            m_camera.drift = {};
            if (m_camera.motion != CameraMotion::Hover) setCameraMotion(CameraMotion::Hover);
            m_settings.flareCamera = false;
            keyframes::applyPose(m_track.path.at(m_track.path.arrival(static_cast<size_t>(ui.keyframeGo))),
                                 m_settings.spin, m_camera.minDistance(m_settings.spin), m_camera);
        }
        if (ui.mergerRequest == 1) m_merger.start(m_camera);
        if (ui.mergerRequest == 2) {
            m_merger.stop();
            m_camera.reset();
        }
        if (ui.burstRequested) {
            const bool dwarf = m_settings.central.kind == central::Kind::WhiteDwarf;
            triggerFlash(dwarf ? 10.0 : 4.0, dwarf ? 300.0 : 30.0, dwarf ? "Nova." : "Type I X-ray burst.");
        }
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
            if (!m_renderer.shaderErrors().empty() || !m_particles.errors().empty() || !m_merger.errors().empty()) {
                std::fprintf(stderr, "%s%s%s\n", m_renderer.shaderErrors().c_str(), m_particles.errors().c_str(),
                             m_merger.errors().c_str());
                exitCode = 2;
            }
            saveScreenshot(options.screenshotPath);
            {
                const central::Object& o = m_settings.central;
                std::printf("  central: %s (%s), M = %.4g Msun, a = %.4f, scene scale %.4g", o.name.c_str(),
                            central::kindName(o.kind), o.massSolar, m_settings.spin, sceneScale());
                if (!o.hasHorizon()) {
                    std::printf(", R = %.6g M = %s, 1 + z = %.6f", o.radius, central::formatLength(o.radius, o.massSolar).c_str(),
                                1.0 + central::surfaceRedshift(o.radius));
                }
                std::printf("\n");
            }
            for (const Body& b : m_bodies.bodies) {
                const auto x = kerr::toCartesian(m_settings.spin, b.state.r, b.state.theta, b.state.phi);
                std::printf("  %-16s %-12s r = %8.3f  theta = %6.3f  E = %.6f  L = %+.5f  Q = %.5f  "
                            "at (%.2f, %.2f, %.2f)\n",
                            b.name.c_str(), kindName(b.kind), b.state.r, b.state.theta, b.constants.E,
                            b.constants.L, b.constants.Q, x[0], x[1], x[2]);
            }
            if (m_particles.active()) std::printf("  debris: %s\n", m_particles.describe().c_str());
            if (m_merger.active()) std::printf("  %s\n", m_merger.describe().c_str());
            if (m_camera.motion != CameraMotion::Hover) {
                glm::vec3 v = m_camera.velocity(m_settings.spin);
                std::printf("  camera: r = %.3f, speed %.3f c, clock rate %.3f, proper time %.2f M%s\n", m_camera.r,
                            glm::length(v), m_camera.clockRate(m_settings.spin), m_camera.properTime,
                            m_camera.ended ? " (the fall has ended)"
                                           : m_camera.insideHorizon(m_settings.spin) ? " (inside the horizon)" : "");
            }
            for (const std::string& e : m_bodies.events) std::printf("  event: %s\n", e.c_str());
            if (!m_ui.cameraStatus.empty()) std::printf("  camera status: %s\n", m_ui.cameraStatus.c_str());
            if (!m_ui.centralStatus.empty()) std::printf("  central status: %s\n", m_ui.centralStatus.c_str());
            std::printf("GPU time per frame: %.2f ms at %dx%d traced\n", m_renderer.gpuMilliseconds(),
                        m_renderer.traceWidth(), m_renderer.traceHeight());
            break;
        }
    }
    shutdown();
    return exitCode;
}
