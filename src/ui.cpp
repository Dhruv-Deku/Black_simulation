#include "ui.h"

#include "blackbody.h"
#include "kerr.h"
#include "renderer.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

constexpr float kDeg = 3.14159265f / 180.0f;
const ImU32 kPlannerColor = IM_COL32(80, 220, 255, 255);
const ImU32 kRayColor = IM_COL32(255, 230, 90, 255);

void help(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

ImU32 toColor(const glm::vec3& c, float alpha = 1.0f) {
    return ImGui::GetColorU32(ImVec4(c.r, c.g, c.b, alpha));
}

// The custom launch's velocity in the ZAMO frame. False if it is not possible here.
bool plannedVelocity(const UiState& ui, double spin, double v[3]) {
    const double r = ui.launchR;
    if (r < kerr::horizonRadius(spin) + 0.2) return false;
    const double vc = r > kerr::photonOrbitRadius(spin, true) * 1.05 ? kerr::circularZamoSpeed(spin, r, true)
                                                                      : 1.0 / std::sqrt(std::max(r - 2.0, 0.5));
    const double speed = std::min(ui.launchSpeed * vc, 0.98);
    const double tilt = ui.launchTiltDeg * kDeg;
    const double dir = ui.launchPrograde ? 1.0 : -1.0;
    v[0] = ui.launchRadial;
    v[1] = -std::sin(tilt) * speed;
    v[2] = dir * std::cos(tilt) * speed;
    const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if (v2 >= 0.995 * 0.995) {
        for (int i = 0; i < 3; ++i) v[i] *= 0.995 / std::sqrt(v2);
    }
    return true;
}

// --- Sections ------------------------------------------------------------------

void statsSection(UiContext& ctx) {
    const Renderer& r = ctx.renderer;
    ImGui::Text("%.0f fps   GPU %.1f ms   trace %dx%d", ctx.fps, r.gpuMilliseconds(), r.traceWidth(),
                r.traceHeight());
    ImGui::TextDisabled("%s", ctx.gpuName);
    ImGui::Text("t = %.1f M", ctx.bodies.time);
    ImGui::SameLine();
    ImGui::TextDisabled("(coordinate time of a distant observer)");
    ImGui::Checkbox("Ray inspector", &ctx.settings.rayInspector);
    help("When on, Ctrl+click or middle-click any pixel to trace the light that reached it: where it "
         "came from, its frequency shift, how often it circled the hole and how long it travelled. "
         "Off, clicks do nothing special and the inspector is hidden. Shortcut: I.");
    if (ctx.settings.rayInspector) {
        ImGui::SameLine();
        ImGui::TextDisabled("Ctrl+click a pixel");
    }
}

void blackHoleSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    ImGui::SliderFloat("Spin a/M", &s.spin, 0.0f, 0.999f, "%.3f");
    help("Angular momentum per unit mass. 0 is a non-rotating (Schwarzschild) hole; values near 1 "
         "are extreme Kerr holes whose rotation drags space around them. Most astrophysical black "
         "holes are thought to spin quickly.");

    ImGui::TextDisabled("Famous black holes:");
    for (int i = 1; i < SceneCount; ++i) {
        if (i > 1) ImGui::SameLine();
        if (ImGui::SmallButton(kSceneNames[i])) ctx.sceneRequested = i;
        ImGui::SetItemTooltip("%s", kSceneDescriptions[i]);
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Reset")) ctx.sceneRequested = SceneDefault;

    const double a = s.spin;
    if (ImGui::BeginTable("radii", 3, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Radius (M)");
        ImGui::TableSetupColumn("Prograde");
        ImGui::TableSetupColumn("Retrograde");
        ImGui::TableHeadersRow();
        auto row = [](const char* name, double pro, double retro, const char* tip) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(name);
            ImGui::SetItemTooltip("%s", tip);
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", pro);
            ImGui::TableNextColumn();
            if (retro > 0) ImGui::Text("%.3f", retro);
        };
        row("Event horizon", kerr::horizonRadius(a), -1, "Nothing, not even light, returns from inside.");
        row("Ergosphere (equator)", kerr::ergosphereRadius(a, 0.5 * kerr::kPi), -1,
            "Inside, frame dragging is so strong that nothing can stay still relative to the stars.");
        row("Photon orbit", kerr::photonOrbitRadius(a, true), kerr::photonOrbitRadius(a, false),
            "Light can circle the hole here (unstably). It forms the bright photon ring.");
        row("ISCO", kerr::iscoRadius(a, true), kerr::iscoRadius(a, false),
            "Innermost stable circular orbit: where the accretion disk ends.");
        ImGui::EndTable();
    }
}

void cameraSection(UiContext& ctx) {
    Camera& cam = ctx.camera;
    Settings& s = ctx.settings;
    UiState& ui = ctx.state;
    const float minR = static_cast<float>(kerr::horizonRadius(s.spin)) + 0.3f;

    // Motion
    const CameraMotion motion = cam.motion;
    if (ImGui::RadioButton("Hover", motion == CameraMotion::Hover)) ctx.cameraMotionRequested = 0;
    ImGui::SetItemTooltip("Held in place: a zero-angular-momentum observer (ZAMO) that rotates with the "
                          "space dragged around by the hole.");
    ImGui::SameLine();
    if (ImGui::RadioButton("Orbit", motion == CameraMotion::Orbit)) ctx.cameraMotionRequested = 1;
    ImGui::SetItemTooltip("Ride a circular orbit at the current distance. Moving near light speed, the "
                          "sky ahead crowds together and turns blue (relativistic aberration).");
    ImGui::SameLine();
    if (ImGui::RadioButton("Free fall", motion == CameraMotion::FreeFall)) ctx.cameraMotionRequested = 2;
    ImGui::SetItemTooltip("Drop from rest straight toward the hole. Watch the shadow swell as you fall, then "
                          "cross the horizon and see the outside universe shrink to a bright patch above you.");

    // Ride along with a body.
    const char* riding = "nobody";
    for (const Body& b : ctx.bodies.bodies) {
        if (motion == CameraMotion::RideAlong && b.slot == cam.rideSlot) riding = b.name.c_str();
    }
    if (ImGui::BeginCombo("Ride along with", riding)) {
        if (ImGui::Selectable("nobody", motion != CameraMotion::RideAlong)) ctx.cameraMotionRequested = 0;
        for (const Body& b : ctx.bodies.bodies) {
            bool selected = motion == CameraMotion::RideAlong && b.slot == cam.rideSlot;
            ImGui::PushID(b.slot);
            if (ImGui::Selectable(b.name.c_str(), selected)) ctx.rideRequested = b.slot;
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    help("Sit on an orbiting body and see what it sees: ride the zoom-whirl star round the hole, or a "
         "plunging star down toward the horizon.");

    if (motion != CameraMotion::Hover) {
        if (cam.insideHorizon(s.spin)) {
            ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.35f, 1.0f), "Inside the event horizon   r = %.3f M", cam.r);
            ImGui::TextDisabled("r now works like time: it can only decrease.");
            ImGui::Text("Your clock: %.2f M", cam.properTime);
        } else {
            glm::vec3 v = cam.velocity(s.spin);
            float speed = glm::length(v);
            ImGui::Text("Speed %.3f c   gamma %.2f   r = %.2f M", speed,
                        1.0f / std::sqrt(std::max(1.0f - speed * speed, 1e-6f)), cam.r);
            ImGui::Text("Your clock: %.1f M (runs at %.0f%% of a distant clock)", cam.properTime,
                        100.0 * cam.clockRate(s.spin));
        }
        if (ctx.slowdown < 0.999) ImGui::TextColored(ImVec4(0.6f, 0.8f, 1.0f, 1.0f), "Slow motion: x%.2f", ctx.slowdown);
        ImGui::Checkbox("Slow motion near the hole", &s.slowNearHole);
        help("Falling in, the last stretch and the whole trip inside the horizon take only a few M; time "
             "slows down automatically so you can watch it.");
        ImGui::Checkbox("Look where you're going", &cam.lookAhead);
        help("Faces the direction of travel. At high speed the whole sky crowds toward the front and "
             "turns blue, while the hole appears off to the side or behind.");
        ImGui::SameLine();
        if (ImGui::SmallButton("Face hole##moving")) {
            cam.lookAhead = false;
            cam.yaw = cam.pitch = 0.0f;
        }
        ImGui::TextDisabled("Drag to look around. Pause (Space) stops the motion.");
    } else {
        ImGui::TextDisabled("Clocks here run at %.0f%% of the rate far away.", 100.0 * cam.clockRate(s.spin));
    }
    if (!ui.cameraStatus.empty()) ImGui::TextWrapped("%s", ui.cameraStatus.c_str());
    if (motion == CameraMotion::Hover) {
        float r = static_cast<float>(cam.r);
        if (ImGui::SliderFloat("Distance", &r, minR, 500.0f, "%.2f M", ImGuiSliderFlags_Logarithmic)) cam.r = r;
        float incl = static_cast<float>(cam.theta) / kDeg;
        if (ImGui::SliderFloat("Inclination", &incl, 1.0f, 179.0f, "%.1f deg")) cam.theta = incl * kDeg;
        help("Angle from the spin axis: 0 looks down on the north pole, 90 is edge-on to the disk.");
        float az = static_cast<float>(std::fmod(cam.phi / kDeg + 720.0, 360.0));
        if (ImGui::SliderFloat("Azimuth", &az, 0.0f, 360.0f, "%.1f deg")) cam.phi = az * kDeg;
    }
    ImGui::SliderFloat("Field of view", &cam.fovDeg, 15.0f, 120.0f, "%.0f deg");
    if (motion == CameraMotion::Hover) {
        ImGui::Checkbox("Auto orbit", &s.autoOrbit);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(-1);
        ImGui::SliderFloat("##orbitspeed", &s.autoOrbitSpeed, -0.5f, 0.5f, "%.2f rad/s");

        auto view = [&](const char* label, double rr, double inclDeg) {
            if (ImGui::Button(label)) {
                cam.r = std::max<double>(rr, minR);
                cam.theta = inclDeg * kDeg;
                cam.yaw = cam.pitch = 0.0f;
            }
        };
        view("Cinematic", 32.0, 84.0);
        ImGui::SameLine();
        view("Edge-on", 30.0, 90.0);
        ImGui::SameLine();
        view("Top-down", 30.0, 2.0);
        ImGui::SameLine();
        view("Close", 9.0, 80.0);
        ImGui::SameLine();
    }
    if (ImGui::Button("Face hole")) cam.yaw = cam.pitch = 0.0f;
}

void diskSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    ImGui::Checkbox("Show disk", &s.diskEnabled);
    ImGui::SameLine();
    if (ImGui::RadioButton("Prograde", s.diskPrograde)) s.diskPrograde = true;
    ImGui::SameLine();
    if (ImGui::RadioButton("Retrograde", !s.diskPrograde)) s.diskPrograde = false;

    ImGui::Checkbox("Inner edge at ISCO", &s.diskInnerAtIsco);
    if (!s.diskInnerAtIsco) {
        ImGui::SliderFloat("Inner radius", &s.diskInner, 1.0f, 15.0f, "%.2f M");
    }
    ImGui::SliderFloat("Outer radius", &s.diskOuter, 4.0f, 60.0f, "%.1f M");
    ImGui::SliderFloat("Peak temperature", &s.diskTemperature, 1500.0f, 40000.0f, "%.0f K", ImGuiSliderFlags_Logarithmic);
    help("Sets the disk's color through the blackbody spectrum. Real disks around stellar black holes "
         "reach millions of kelvin; here it is kept in the visible range so the colors mean something.");
    ImGui::SliderFloat("Brightness", &s.diskBrightness, 0.05f, 10.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
    ImGui::SliderFloat("Opacity", &s.diskOpacity, 0.0f, 1.0f, "%.2f");
    ImGui::SliderFloat("Turbulence", &s.diskTurbulence, 0.0f, 1.0f, "%.2f");
    ImGui::Combo("Relativistic light", &s.relativity, "Full (Doppler + gravity)\0Beaming only\0Off (Interstellar style)\0");
    help("Full: approaching gas looks hotter and brighter, receding gas cooler and dimmer, and light "
         "climbing out of the hole's gravity is redshifted.\nBeaming only: brightness changes but colors "
         "do not.\nOff: no frequency shifts, like the film Interstellar, which left them out so audiences "
         "would not be confused.");

    ImGui::Checkbox("Relativistic jets", &s.jetsEnabled);
    help("Twin beams of plasma fired along the spin axis at nearly the speed of light, as from M87. The "
         "jet coming toward you is boosted by Doppler beaming; the one going away nearly vanishes. "
         "Bright knots travel outward along them.");
    if (s.jetsEnabled) {
        ImGui::SliderFloat("Jet power", &s.jetPower, 0.05f, 10.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Jet speed", &s.jetSpeed, 0.0f, 0.99f, "%.2f c");
        ImGui::SliderFloat("Jet opening", &s.jetOpeningDeg, 1.0f, 20.0f, "%.1f deg");
    }
}

void skySection(UiContext& ctx) {
    Settings& s = ctx.settings;
    UiState& ui = ctx.state;
    ImGui::Combo("Background", &s.skyMode, "Stars and galaxy\0Lensing grid\0Image\0Black\0");
    help("The lensing grid makes distortion easy to read: every square is 15 degrees, and the color "
         "shows which quarter of the sky it came from. Look for the multiple images near the shadow.");
    if (s.skyMode == 0) {
        ImGui::SliderFloat("Stars", &s.starBrightness, 0.0f, 4.0f, "%.2f");
        ImGui::SliderFloat("Milky Way", &s.galaxyBrightness, 0.0f, 4.0f, "%.2f");
    }
    if (s.skyMode == 2) {
        ImGui::InputTextWithHint("##skypath", "equirectangular .jpg/.png/.hdr", ui.skyPath, sizeof(ui.skyPath));
        ImGui::SameLine();
        if (ImGui::Button("Load")) {
            std::string error;
            ui.skyStatus = ctx.renderer.loadSkyImage(ui.skyPath, error) ? "Loaded." : "Failed: " + error;
        }
        if (!ctx.renderer.hasSkyImage()) ImGui::TextDisabled("No image loaded: showing stars instead.");
        if (!ui.skyStatus.empty()) ImGui::TextWrapped("%s", ui.skyStatus.c_str());
        ImGui::SliderFloat("Image brightness", &s.skyImageBrightness, 0.0f, 4.0f, "%.2f");
    }
}

void presetButtons(UiContext& ctx, bool objects) {
    const float buttonWidth = (ImGui::GetContentRegionAvail().x - 2.0f * ImGui::GetStyle().ItemSpacing.x) / 3.0f;
    int column = 0;
    for (const PresetInfo& p : presetList()) {
        if (p.object != objects) continue;
        if (column++ % 3 != 0) ImGui::SameLine();
        if (ImGui::Button(p.label, ImVec2(buttonWidth, 0))) ctx.bodies.addPreset(ctx.settings.spin, p.preset);
        ImGui::SetItemTooltip("%s", p.description);
    }
    if (objects) {
        if (column++ % 3 != 0) ImGui::SameLine();
        if (ImGui::Button("Einstein ring", ImVec2(buttonWidth, 0))) ctx.einsteinRingRequested = true;
        ImGui::SetItemTooltip("Puts a star exactly behind the hole, as seen from the camera. Its light "
                              "comes around every side at once and forms a ring.");
    }
}

// --- Launch planner ----------------------------------------------------------------

// Recomputes the planner's prediction when its inputs change.
void updatePlanner(UiContext& ctx) {
    UiState& ui = ctx.state;
    const double a = ctx.settings.spin;
    char key[256];
    std::snprintf(key, sizeof(key), "%.6f %.4f %.3f %.3f %.4f %.3f %.4f %d", a, ui.launchR, ui.launchPolarDeg,
                  ui.launchAzimuthDeg, ui.launchSpeed, ui.launchTiltDeg, ui.launchRadial, ui.launchPrograde ? 1 : 0);
    if (ui.plannerKey == key) return;
    ui.plannerKey = key;
    ui.plannerValid = plannedLaunch(ui, a, ui.plannerConstants, ui.plannerState);
    if (!ui.plannerValid) return;
    ui.plannerSpin = a;
    ui.plannerAnalysis = trajectory::classify(a, ui.plannerConstants, ui.plannerState);
    ui.plannerPath = trajectory::predict(a, ui.plannerConstants, ui.plannerState, 3000.0, 0.25, 8000);
}

// The effective potential V(r): the body can only be where its energy E is above it.
void potentialPlot(const UiState& ui, double a) {
    const kerr::Constants& c = ui.plannerConstants;
    const trajectory::Analysis& an = ui.plannerAnalysis;
    const double rh = kerr::horizonRadius(a);
    const double r0 = ui.plannerState.r;
    double rMax = std::max(40.0, 1.6 * r0);
    if (an.hasOuter) rMax = std::max(rMax, 1.3 * an.rOuter);
    rMax = std::min(rMax, 2000.0);
    const double rMin = rh * 1.001;

    const int n = 240;
    std::vector<double> rs(n), vs(n);
    double vMin = 1e30, vPeak = -1e30;
    for (int i = 0; i < n; ++i) {
        rs[i] = rMin * std::pow(rMax / rMin, static_cast<double>(i) / (n - 1));
        vs[i] = trajectory::effectivePotential(a, c.L, c.Q, rs[i]);
        vMin = std::min(vMin, vs[i]);
        vPeak = std::max(vPeak, vs[i]);
    }
    const double E = c.E;
    double yHi = std::max(E, std::min(vPeak, E + 0.25)) + 0.03;
    double yLo = std::min(vMin, E);
    yLo -= 0.25 * (yHi - yLo);  // room below the well for the energy line

    const float em = ImGui::GetFontSize();
    ImVec2 origin = ImGui::GetCursorScreenPos();
    ImVec2 size(ImGui::GetContentRegionAvail().x, 8.5f * em);
    ImGui::InvisibleButton("potential", size);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 corner(origin.x + size.x, origin.y + size.y);
    dl->AddRectFilled(origin, corner, IM_COL32(10, 11, 16, 255));
    dl->PushClipRect(origin, corner, true);
    auto px = [&](double r) {
        return origin.x + static_cast<float>(std::log(r / rMin) / std::log(rMax / rMin)) * size.x;
    };
    auto py = [&](double v) { return corner.y - static_cast<float>((v - yLo) / (yHi - yLo)) * size.y; };

    // Forbidden regions (V > E) shaded.
    for (int i = 1; i < n; ++i) {
        if (vs[i] > E || vs[i - 1] > E) {
            dl->AddRectFilled(ImVec2(px(rs[i - 1]), origin.y), ImVec2(px(rs[i]) + 1.0f, corner.y), IM_COL32(90, 25, 30, 120));
        }
    }
    // ISCO and photon orbit markers.
    const double isco = kerr::iscoRadius(a, true);
    if (isco < rMax) {
        dl->AddLine(ImVec2(px(isco), origin.y), ImVec2(px(isco), corner.y), IM_COL32(80, 220, 120, 70));
    }
    // V(r) and E.
    std::vector<ImVec2> curve(n);
    for (int i = 0; i < n; ++i) curve[i] = ImVec2(px(rs[i]), py(vs[i]));
    dl->AddPolyline(curve.data(), n, IM_COL32(230, 230, 240, 255), 0, 1.5f);
    dl->AddLine(ImVec2(origin.x, py(E)), ImVec2(corner.x, py(E)), kPlannerColor, 1.5f);
    // Turning points and the start.
    if (an.hasInner) dl->AddCircleFilled(ImVec2(px(an.rInner), py(E)), 3.5f, IM_COL32(255, 255, 255, 255));
    if (an.hasOuter && an.rOuter < rMax) dl->AddCircleFilled(ImVec2(px(an.rOuter), py(E)), 3.5f, IM_COL32(255, 255, 255, 255));
    dl->AddCircleFilled(ImVec2(px(r0), py(E)), 5.0f, kPlannerColor);
    dl->PopClipRect();
    dl->AddText(ImVec2(origin.x + 4, origin.y + 2), IM_COL32(170, 170, 180, 255), "V(r)");
    dl->AddText(ImVec2(corner.x - 5.0f * em, py(E) + 2.0f), kPlannerColor, "energy E");
    char label[64];
    std::snprintf(label, sizeof(label), "horizon ... r (log) ... %.0f M", rMax);
    const float labelWidth = ImGui::CalcTextSize(label).x;
    dl->AddText(ImVec2(corner.x - labelWidth - 4, origin.y + 2), IM_COL32(150, 150, 160, 255), label);
}

void plannerSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    BodySystem& sys = ctx.bodies;
    UiState& ui = ctx.state;
    const double a = s.spin;
    ui.plannerOpen = true;

    ImGui::TextWrapped("Set a start and a velocity, or drag on the orbit map (top view) from the start "
                       "point in the direction to throw. The verdict and the dashed path update live.");
    ImGui::Combo("Type", &ui.launchKind, "Star\0Rocky planet\0Gas giant\0Pulsar\0Star (tidally disrupted)\0");
    ImGui::SliderFloat("Start radius", &ui.launchR, 2.0f, 80.0f, "%.2f M", ImGuiSliderFlags_Logarithmic);
    ImGui::SliderFloat("Polar angle", &ui.launchPolarDeg, 1.0f, 179.0f, "%.0f deg");
    ImGui::SliderFloat("Azimuth##launch", &ui.launchAzimuthDeg, 0.0f, 360.0f, "%.0f deg");
    ImGui::SliderFloat("Speed", &ui.launchSpeed, 0.0f, 2.0f, "%.3f x circular");
    help("Sideways speed as a multiple of the circular-orbit speed at the start radius. Below 1 the "
         "orbit dips inward; around 1.4 it escapes.");
    ImGui::SliderFloat("Tilt", &ui.launchTiltDeg, -90.0f, 90.0f, "%.0f deg");
    help("Tilts the velocity out of the equatorial direction. 90 sends the body over the poles.");
    ImGui::SliderFloat("Radial speed", &ui.launchRadial, -0.95f, 0.95f, "%.3f c");
    ImGui::Checkbox("Prograde", &ui.launchPrograde);
    const BodyKind kind = static_cast<BodyKind>(ui.launchKind);
    const bool planet = kind == BodyKind::RockyPlanet || kind == BodyKind::GasGiant;
    if (!planet) {
        ImGui::SliderFloat("Temperature", &ui.launchTemperature, 2000.0f, 40000.0f, "%.0f K", ImGuiSliderFlags_Logarithmic);
    }
    ImGui::SliderFloat("Size", &ui.launchRadius, 0.1f, 2.0f, "%.2f M");
    if (kind != BodyKind::TidalStar) {
        ImGui::SliderFloat("Mass", &ui.launchMass, 1e-7f, 0.05f, "%.2g M", ImGuiSliderFlags_Logarithmic);
        help("In units of the hole's mass. A Sun-like star around a million-Sun hole is 1e-6: its pull "
             "on anything else is negligible. Raise it (to ~0.01) and, with gravity between bodies on, it "
             "can hold moons and capture or fling other bodies.");
    }

    // Verdict
    updatePlanner(ctx);
    if (ui.plannerValid) {
        const trajectory::Analysis& an = ui.plannerAnalysis;
        const trajectory::Path& path = ui.plannerPath;
        ImVec4 color(0.4f, 0.9f, 0.5f, 1.0f);
        char verdict[256];
        switch (an.fate) {
            case trajectory::Fate::Plunges:
                color = ImVec4(1.0f, 0.4f, 0.35f, 1.0f);
                if (path.fate == trajectory::Fate::Plunges) {
                    std::snprintf(verdict, sizeof(verdict), "Falls into the black hole after %.0f M.", path.duration);
                } else {
                    std::snprintf(verdict, sizeof(verdict), "Falls into the black hole%s.",
                                  an.hasOuter ? " after swinging out to r = " : "");
                    if (an.hasOuter) {
                        std::snprintf(verdict, sizeof(verdict), "Swings out to r = %.1f M, then falls in.", an.rOuter);
                    }
                }
                break;
            case trajectory::Fate::Escapes:
                color = ImVec4(1.0f, 0.85f, 0.3f, 1.0f);
                std::snprintf(verdict, sizeof(verdict), "Escapes to infinity%s (E = %.4f).",
                              an.hasInner ? " after swinging past the hole" : "", ui.plannerConstants.E);
                break;
            case trajectory::Fate::Circular:
                if (!an.stableCircular) color = ImVec4(1.0f, 0.65f, 0.3f, 1.0f);
                std::snprintf(verdict, sizeof(verdict), "%s", an.stableCircular
                                  ? "Stable circular orbit."
                                  : "Circular, but unstable: the slightest nudge sends it in or out.");
                break;
            case trajectory::Fate::Bound:
                std::snprintf(verdict, sizeof(verdict), "Bound orbit between r = %.2f and %.2f M.", an.rInner,
                              an.rOuter);
                break;
        }
        ImGui::PushStyleColor(ImGuiCol_Text, color);
        ImGui::TextWrapped("%s", verdict);
        ImGui::PopStyleColor();
        ImGui::TextDisabled("E = %.5f   L = %+.4f   Q = %.4f   tilt %.1f deg", ui.plannerConstants.E,
                            ui.plannerConstants.L, ui.plannerConstants.Q, an.inclinationDeg);
        if (an.fate == trajectory::Fate::Bound && path.periapses >= 2) {
            ImGui::TextDisabled("The closest point advances %.1f deg every orbit.", path.precessionDeg);
        }
        if (ui.plannerJointValid) {
            const trajectory::Path& joint = ui.plannerJointPath;
            char among[256];
            if (joint.fate == trajectory::Fate::Plunges) {
                std::snprintf(among, sizeof(among), "With the other bodies' pull: falls in after %.0f M.", joint.duration);
            } else if (joint.fate == trajectory::Fate::Escapes) {
                std::snprintf(among, sizeof(among), "With the other bodies' pull: escapes (gone after %.0f M).",
                              joint.duration);
            } else {
                std::snprintf(among, sizeof(among),
                              "With the other bodies' pull: still in orbit after %.0f M, between r = %.1f and %.1f M.",
                              joint.duration, joint.rMin, joint.rMax);
            }
            ImGui::TextWrapped("%s", among);
            help("The verdict above is exact for a body that only feels the hole. This line comes from "
                 "simulating it together with the other bodies (gravity between bodies is on), which can "
                 "change its fate. The dashed path shows this simulation.");
        }
        potentialPlot(ui, a);
        help("Effective potential: the least energy a body with this angular momentum needs to be at each "
             "radius. It can only go where its energy line (blue) is above the curve; the dots are the turning "
             "points. A hump the line clears means nothing stops it falling in.");
    } else {
        ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.35f, 1.0f), "Cannot start there: too close to the horizon.");
    }

    auto finish = [&]() {
        Body& b = sys.bodies.back();
        b.temperature = planet ? 0.0f : ui.launchTemperature;
        b.radius = ui.launchRadius;
        if (kind == BodyKind::RockyPlanet) b.spinRate = 0.06f;
        if (kind == BodyKind::GasGiant) b.spinRate = 0.1f;
        if (kind == BodyKind::Pulsar) b.spinRate = 2.0f;
        b.mass = ui.launchMass;
        if (kind == BodyKind::TidalStar) {
            b.tidalRadius = static_cast<float>(std::max(15.0, 0.6 * ui.launchR));
            b.mass = std::pow(b.radius / b.tidalRadius, 3.0);  // r_t = R (M / m)^(1/3)
        }
    };
    const bool enter = !ImGui::GetIO().WantTextInput && ImGui::IsKeyPressed(ImGuiKey_Enter, false);
    if (ImGui::Button("Launch (Enter)") || enter) {
        double v[3];
        if (plannedVelocity(ui, a, v) &&
            sys.addFromZamo(a, ui.launchR, ui.launchPolarDeg * kDeg, ui.launchAzimuthDeg * kDeg, v, kind)) {
            finish();
        }
    }
    if (kind == BodyKind::TidalStar) {
        ImGui::SameLine();
        ImGui::TextDisabled("Tidal radius %.0f M: torn apart well inside it.", std::max(15.0, 0.6 * ui.launchR));
    }
    ImGui::SeparatorText("Or by turning points (equatorial)");
    ImGui::SliderFloat("Periapsis", &ui.launchRp, 1.0f, 40.0f, "%.2f M");
    ImGui::SliderFloat("Apoapsis", &ui.launchRa, 2.0f, 80.0f, "%.2f M");
    if (ImGui::Button("Launch bound orbit")) {
        if (sys.addEquatorial(a, ui.launchRp, ui.launchRa, ui.launchPrograde, ui.launchAzimuthDeg * kDeg, kind)) {
            finish();
        }
    }
}

void objectsSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    BodySystem& sys = ctx.bodies;

    ImGui::Checkbox("Paused", &s.paused);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-1);
    ImGui::SliderFloat("##timescale", &s.timeScale, 0.1f, 500.0f, "%.1f M per second", ImGuiSliderFlags_Logarithmic);
    ImGui::Checkbox("Light travel time", &s.lightDelay);
    help("Light takes time to reach you, so everything is seen where it was when its light set out: "
         "fast stars lag behind and stretch, lensed images of the same star show it at different "
         "moments, and a pulsar's lensed images flash out of step.");
    ImGui::Checkbox("Gravity between bodies", &s.bodyGravity);
    help("Bodies pull on each other (and on debris gas), each feeling the others' Newtonian pull in its "
         "own rest frame on top of the hole's exact gravity. Most bodies are far too light to matter; the "
         "binary star and planet-with-moon presets are heavy enough to show it.");
    ImGui::Checkbox("Paths in view", &s.showPaths);
    help("Draws each body's trail, and its predicted future (dashed), into the scene. They are ray "
         "traced like everything else, so a path behind the hole appears bent around it.");
    if (s.showPaths) {
        ImGui::SameLine();
        ImGui::Checkbox("Predicted", &s.showPredictions);
        if (s.showPredictions) {
            ImGui::SliderFloat("Predict ahead", &s.predictionTime, 50.0f, 3000.0f, "%.0f M", ImGuiSliderFlags_Logarithmic);
        }
        ImGui::SliderFloat("Line width", &s.pathWidth, 0.5f, 4.0f, "%.1f px");
    }

    ImGui::SeparatorText("Objects");
    presetButtons(ctx, true);
    ImGui::SeparatorText("Orbits");
    presetButtons(ctx, false);

    ctx.state.plannerOpen = false;
    if (ImGui::TreeNodeEx("Launch planner", ImGuiTreeNodeFlags_DefaultOpen)) {
        plannerSection(ctx);
        ImGui::TreePop();
    }

    if (ctx.particles.active()) {
        ImGui::SeparatorText("Debris");
        ImGui::Text("%d particles", ctx.particles.alive());
        ImGui::SameLine();
        if (ImGui::SmallButton("Clear")) ctx.particles.clear();
        ImGui::SliderFloat("Glow", &s.debrisBrightness, 0.1f, 20.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Thickness", &s.debrisOpacity, 0.005f, 3.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
        help("How strongly the gas absorbs light: low is a faint translucent glow, high is dense and opaque.");
        ImGui::Checkbox("Gas physics", &s.gasPhysics);
        help("Pressure pushes the gas apart; where streams of gas cross they collide and shock, turning "
             "their motion into heat (this is what slowly turns tidal debris into a disk); squeezed gas "
             "heats up, expanding gas cools, and hot gas glows its heat away. Off, every parcel just "
             "follows its own orbit.");
        if (s.gasPhysics) {
            ImGui::SliderFloat("Shocks", &s.gasCollisions, 0.0f, 3.0f, "%.2f");
            help("How readily gas parcels in the same place collide: 0 lets crossing streams pass through each "
                 "other.");
            ImGui::SliderFloat("Cooling time", &s.gasCoolingTime, 5.0f, 1000.0f, "%.0f M", ImGuiSliderFlags_Logarithmic);
            help("How long gas at 10,000 K takes to radiate its heat away; hotter gas cools faster (~T^-3).");
            ImGui::SliderFloat("Coolest", &s.gasFloorTemperature, 1000.0f, 15000.0f, "%.0f K", ImGuiSliderFlags_Logarithmic);
            help("The light around the hole keeps the gas from cooling below this.");
        }
    }

    if (!sys.bodies.empty()) {
        ImGui::SeparatorText("Bodies");
        int removeIndex = -1;
        if (ImGui::BeginTable("bodies", 7, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Name");
            ImGui::TableSetupColumn("r");
            ImGui::TableSetupColumn("E");
            ImGui::TableSetupColumn("L");
            ImGui::TableSetupColumn("Q");
            ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed);
            ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed);
            ImGui::TableHeadersRow();
            for (int i = 0; i < static_cast<int>(sys.bodies.size()); ++i) {
                const Body& b = sys.bodies[i];
                ImGui::PushID(i);
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                glm::vec3 hue = displayColor(b);
                ImGui::ColorButton("##hue", ImVec4(hue.r, hue.g, hue.b, 1.0f), ImGuiColorEditFlags_NoTooltip,
                                   ImVec2(ImGui::GetFontSize() * 0.8f, ImGui::GetFontSize() * 0.8f));
                ImGui::SameLine();
                ImGui::TextUnformatted(b.name.c_str());
                ImGui::SetItemTooltip("%s, mass %.2g M%s", kindName(b.kind), b.mass,
                                      b.inside ? " (inside the horizon)" : "");
                ImGui::TableNextColumn();
                ImGui::Text("%.2f", b.state.r);
                ImGui::TableNextColumn();
                ImGui::Text("%.4f", b.constants.E);
                ImGui::SetItemTooltip("Energy per unit rest mass, conserved. Below 1 means bound.");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", b.constants.L);
                ImGui::SetItemTooltip("Angular momentum about the spin axis, conserved.");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", b.constants.Q);
                ImGui::SetItemTooltip("Carter constant, conserved. Zero for equatorial orbits.");
                ImGui::TableNextColumn();
                if (ImGui::SmallButton("ride")) ctx.rideRequested = b.slot;
                ImGui::SetItemTooltip("Put the camera on this body.");
                ImGui::TableNextColumn();
                if (ImGui::SmallButton("x")) removeIndex = i;
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        if (removeIndex >= 0) sys.remove(static_cast<size_t>(removeIndex));
        if (ImGui::Button("Remove all")) {
            while (!sys.bodies.empty()) sys.remove(sys.bodies.size() - 1);
        }
    }
    for (const std::string& e : sys.events) ImGui::TextDisabled("%s", e.c_str());
}

void videoSection(UiContext& ctx) {
    UiState& ui = ctx.state;
    if (ctx.recording) {
        float progress = ctx.recordTotal > 0 ? static_cast<float>(ctx.recordFrame) / ctx.recordTotal : 0.0f;
        char label[64];
        std::snprintf(label, sizeof(label), "frame %d / %d", ctx.recordFrame, ctx.recordTotal);
        ImGui::ProgressBar(progress, ImVec2(-1, 0), label);
        if (ImGui::Button("Stop (Esc)")) ctx.stopRecordingRequested = true;
    } else {
        const char* sizes[] = {kVideoSizes[0].label, kVideoSizes[1].label, kVideoSizes[2].label, kVideoSizes[3].label};
        ImGui::Combo("Size", &ui.videoSize, sizes, 4);
        ImGui::Combo("Frame rate", &ui.videoFps, "30 fps\0" "60 fps\0");
        ImGui::SliderFloat("Length", &ui.videoSeconds, 1.0f, 120.0f, "%.0f s", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Quality", &ui.videoQuality, 0.5f, 2.0f, "%.2f x");
        help("Ray tracing resolution relative to the video size. Above 1 smooths edges (supersampling) "
             "and takes longer. Recording runs slower than real time if needed; the video still plays "
             "smoothly.");
        ImGui::Combo("Camera", &ui.videoPath, "As you fly it\0One orbit around the hole\0Slow approach\0");
        if (ImGui::Button("Record video")) ctx.recordRequested = true;
        ImGui::SameLine();
        ImGui::TextDisabled("saved in videos\\");
    }
    if (!ui.videoStatus.empty()) ImGui::TextWrapped("%s", ui.videoStatus.c_str());
}

void renderSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    ImGui::SliderFloat("Resolution", &s.renderScale, 0.25f, 2.0f, "%.2f x");
    help("Fraction of the window resolution that is ray traced. Lower is faster; above 1 supersamples.");
    ImGui::SliderFloat("Step size", &s.stepSize, 0.01f, 0.15f, "%.3f", ImGuiSliderFlags_Logarithmic);
    help("Integration accuracy: smaller steps trace the strongly bent rays near the photon ring more "
         "precisely, at a higher cost.");
    ImGui::SliderInt("Max steps", &s.maxSteps, 100, 3000);
    ImGui::SliderFloat("Exposure", &s.exposure, 0.05f, 10.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
    ImGui::SliderFloat("Bloom", &s.bloomStrength, 0.0f, 0.3f, "%.3f");
    ImGui::Combo("Tone curve", &s.tonemapper, "ACES filmic\0Reinhard\0Clamp\0");
    ImGui::Checkbox("VSync", &s.vsync);
    ImGui::SameLine();
    bool heat = s.debugView == 1;
    if (ImGui::Checkbox("Step heat map", &heat)) s.debugView = heat ? 1 : 0;
    help("Shows how many integration steps each pixel needed. The photon ring is the most expensive part.");
    if (ImGui::Button("Save screenshot (F12)")) ctx.screenshotRequested = true;
    if (!ctx.state.screenshotStatus.empty()) {
        ImGui::SameLine();
        ImGui::TextDisabled("%s", ctx.state.screenshotStatus.c_str());
    }
}

// --- Orbit map -------------------------------------------------------------------

void orbitMap(UiContext& ctx) {
    UiState& ui = ctx.state;
    const Settings& s = ctx.settings;
    const double a = s.spin;

    ImGui::SetNextWindowSize(ImVec2(360, 400), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 380, 20), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Orbit map", &ui.showOrbitMap)) {
        ImGui::End();
        return;
    }
    if (ImGui::RadioButton("From above", ui.mapView == 0)) ui.mapView = 0;
    ImGui::SameLine();
    if (ImGui::RadioButton("From the side", ui.mapView == 1)) ui.mapView = 1;
    ImGui::SameLine();
    ImGui::TextDisabled("(scroll to zoom)");

    ImVec2 origin = ImGui::GetCursorScreenPos();
    ImVec2 size = ImGui::GetContentRegionAvail();
    size.x = std::max(size.x, 50.0f);
    size.y = std::max(size.y, 50.0f);
    ImGui::InvisibleButton("map", size);
    const bool hovered = ImGui::IsItemHovered();
    if (hovered && ImGui::GetIO().MouseWheel != 0.0f) {
        ui.mapExtent = std::clamp(ui.mapExtent * std::pow(0.85f, ImGui::GetIO().MouseWheel), 3.0f, 400.0f);
    }

    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 maxCorner(origin.x + size.x, origin.y + size.y);
    dl->AddRectFilled(origin, maxCorner, IM_COL32(8, 9, 14, 255));
    dl->PushClipRect(origin, maxCorner, true);

    const ImVec2 center(origin.x + size.x * 0.5f, origin.y + size.y * 0.5f);
    const float scale = 0.5f * std::min(size.x, size.y) / ui.mapExtent;
    const bool side = ui.mapView == 1;
    // Project a Cartesian point (z = spin axis) onto the map.
    auto project = [&](float x, float y, float z) {
        float u = x, v = side ? z : y;
        return ImVec2(center.x + u * scale, center.y - v * scale);
    };
    auto projectPoint = [&](const trajectory::Point& p) {
        return project(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
    };
    auto dashedPath = [&](const std::vector<trajectory::Point>& pts, ImU32 color, float thickness) {
        float along = 0.0f;
        for (size_t i = 1; i < pts.size(); ++i) {
            ImVec2 q0 = projectPoint(pts[i - 1]), q1 = projectPoint(pts[i]);
            along += std::hypot(q1.x - q0.x, q1.y - q0.y);
            if (std::fmod(along, 10.0f) < 6.0f) dl->AddLine(q0, q1, color, thickness);
        }
    };

    // Aiming a launch: press on the top view to set the start, drag to set the velocity.
    if (ui.plannerOpen && !side) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const float mx = (mouse.x - center.x) / scale, my = -(mouse.y - center.y) / scale;
        if (ImGui::IsItemActivated()) {
            auto [r, theta, phi] = kerr::fromCartesian(a, mx, my, 0.0);
            if (r > kerr::horizonRadius(a) + 0.3) {
                ui.mapAiming = true;
                ui.launchR = static_cast<float>(r);
                ui.launchPolarDeg = 90.0f;
                ui.launchAzimuthDeg = static_cast<float>(std::fmod(phi / kDeg + 360.0, 360.0));
                ui.launchTiltDeg = 0.0f;
            }
        }
        if (ui.mapAiming && ImGui::IsItemActive()) {
            const double phi = ui.launchAzimuthDeg * kDeg;
            const double sx = ui.launchR * std::cos(phi), sy = ui.launchR * std::sin(phi);
            // 0.03 c of speed per M of drag.
            double vx = (mx - sx) * 0.03, vy = (my - sy) * 0.03;
            double speed = std::hypot(vx, vy);
            if (speed > 0.95) {
                vx *= 0.95 / speed;
                vy *= 0.95 / speed;
            }
            const double vr = vx * std::cos(phi) + vy * std::sin(phi);
            const double vphi = -vx * std::sin(phi) + vy * std::cos(phi);
            const double r = ui.launchR;
            const double vc = r > kerr::photonOrbitRadius(a, true) * 1.05 ? kerr::circularZamoSpeed(a, r, true)
                                                                          : 1.0 / std::sqrt(std::max(r - 2.0, 0.5));
            ui.launchRadial = static_cast<float>(vr);
            ui.launchSpeed = static_cast<float>(std::abs(vphi) / vc);
            ui.launchPrograde = vphi >= 0.0;
        }
        if (!ImGui::IsItemActive()) ui.mapAiming = false;
    }

    // Accretion disk
    if (s.diskEnabled) {
        double inner = s.diskInnerAtIsco ? kerr::iscoRadius(a, s.diskPrograde) : s.diskInner;
        ImU32 diskCol = IM_COL32(255, 150, 60, 50);
        if (side) {
            float x0 = static_cast<float>(std::sqrt(inner * inner + a * a));
            float x1 = static_cast<float>(std::sqrt(s.diskOuter * s.diskOuter + a * a));
            dl->AddLine(project(x0, 0, 0), project(x1, 0, 0), IM_COL32(255, 150, 60, 140), 3.0f);
            dl->AddLine(project(-x0, 0, 0), project(-x1, 0, 0), IM_COL32(255, 150, 60, 140), 3.0f);
        } else {
            float rIn = static_cast<float>(std::sqrt(inner * inner + a * a)) * scale;
            float rOut = static_cast<float>(std::sqrt(s.diskOuter * s.diskOuter + a * a)) * scale;
            dl->AddCircle(center, 0.5f * (rIn + rOut), diskCol, 128, rOut - rIn);
        }
    }

    // Closed curve of constant Boyer-Lindquist radius r(theta), seen from the chosen side.
    auto shell = [&](auto radiusAt, ImU32 col, float thickness, bool filled) {
        std::vector<ImVec2> pts;
        const int n = 96;
        for (int i = 0; i < n; ++i) {
            double t = 2.0 * kerr::kPi * i / n;
            if (side) {
                double r = radiusAt(std::abs(std::remainder(t, 2.0 * kerr::kPi)));
                auto p = kerr::toCartesian(a, r, t, 0.0);
                pts.push_back(project(static_cast<float>(p[0]), 0.0f, static_cast<float>(p[2])));
            } else {
                double r = radiusAt(0.5 * kerr::kPi);
                auto p = kerr::toCartesian(a, r, 0.5 * kerr::kPi, t);
                pts.push_back(project(static_cast<float>(p[0]), static_cast<float>(p[1]), 0.0f));
            }
        }
        if (filled) dl->AddConvexPolyFilled(pts.data(), n, col);
        else dl->AddPolyline(pts.data(), n, col, ImDrawFlags_Closed, thickness);
    };
    const double rh = kerr::horizonRadius(a);
    shell([&](double th) { return kerr::ergosphereRadius(a, th); }, IM_COL32(120, 90, 200, 200), 1.5f, false);
    shell([&](double) { return rh; }, IM_COL32(0, 0, 0, 255), 0.0f, true);
    shell([&](double) { return rh; }, IM_COL32(200, 200, 220, 255), 1.5f, false);
    if (!side) {
        float isco = static_cast<float>(std::sqrt(std::pow(kerr::iscoRadius(a, true), 2) + a * a)) * scale;
        float iscoR = static_cast<float>(std::sqrt(std::pow(kerr::iscoRadius(a, false), 2) + a * a)) * scale;
        float photon = static_cast<float>(std::sqrt(std::pow(kerr::photonOrbitRadius(a, true), 2) + a * a)) * scale;
        dl->AddCircle(center, isco, IM_COL32(80, 220, 120, 160), 96, 1.0f);
        dl->AddCircle(center, iscoR, IM_COL32(220, 140, 80, 110), 96, 1.0f);
        dl->AddCircle(center, photon, IM_COL32(255, 255, 140, 120), 96, 1.0f);
        // Rotation direction of the hole.
        float rr = static_cast<float>(rh) * scale + 10.0f;
        dl->PathArcTo(center, rr, -0.5f, 0.9f, 16);
        dl->PathStroke(IM_COL32(200, 200, 220, 120), 0, 1.5f);
    } else {
        float tip = static_cast<float>(rh) * scale + 16.0f;
        dl->AddLine(ImVec2(center.x, center.y - tip + 12.0f), ImVec2(center.x, center.y - tip - 6.0f),
                    IM_COL32(200, 200, 220, 160), 1.5f);
        dl->AddTriangleFilled(ImVec2(center.x - 4, center.y - tip - 4), ImVec2(center.x + 4, center.y - tip - 4),
                              ImVec2(center.x, center.y - tip - 11), IM_COL32(200, 200, 220, 160));
    }

    // Bodies, their trails and predicted paths; the body nearest the mouse is labeled.
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    const Body* nearest = nullptr;
    ImVec2 nearestPos;
    float nearestDist2 = 20.0f * 20.0f;
    for (size_t k = 0; k < ctx.bodies.bodies.size(); ++k) {
        const Body& b = ctx.bodies.bodies[k];
        glm::vec3 hue = displayColor(b);
        const size_t n = b.trail.size();
        for (size_t i = 1; i < n; ++i) {
            float age = static_cast<float>(i) / static_cast<float>(n);
            const glm::vec3& p0 = b.trail[i - 1];
            const glm::vec3& p1 = b.trail[i];
            dl->AddLine(project(p0.x, p0.y, p0.z), project(p1.x, p1.y, p1.z), toColor(hue, 0.08f + 0.6f * age), 1.2f);
        }
        if (s.showPredictions && ctx.predictions && b.slot >= 0 &&
            static_cast<size_t>(b.slot) < ctx.predictions->size()) {
            dashedPath((*ctx.predictions)[b.slot].points, toColor(hue, 0.45f), 1.0f);
        }
        auto p = kerr::toCartesian(a, b.state.r, b.state.theta, b.state.phi);
        ImVec2 q = project(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
        dl->AddCircleFilled(q, std::max(3.0f, b.radius * scale), toColor(hue));
        float d2 = (q.x - mouse.x) * (q.x - mouse.x) + (q.y - mouse.y) * (q.y - mouse.y);
        if (hovered && d2 < nearestDist2) {
            nearestDist2 = d2;
            nearest = &b;
            nearestPos = q;
        }
    }

    // The planned launch: its path, start point and velocity arrow.
    if (ui.plannerOpen && ui.plannerValid) {
        dashedPath(ui.plannerJointValid ? ui.plannerJointPath.points : ui.plannerPath.points, kPlannerColor, 1.6f);
        const kerr::State& st = ui.plannerState;
        auto p0 = kerr::toCartesian(a, st.r, st.theta, st.phi);
        ImVec2 q0 = projectPoint(p0);
        dl->AddCircle(q0, 5.0f, kPlannerColor, 12, 2.0f);
        if (ui.plannerPath.points.size() > 3) {
            ImVec2 q1 = projectPoint(ui.plannerPath.points[std::min<size_t>(3, ui.plannerPath.points.size() - 1)]);
            float dx = q1.x - q0.x, dy = q1.y - q0.y, len = std::max(std::hypot(dx, dy), 1e-3f);
            ImVec2 tip(q0.x + dx / len * 22.0f, q0.y + dy / len * 22.0f);
            dl->AddLine(q0, tip, kPlannerColor, 2.0f);
            dl->AddTriangleFilled(ImVec2(tip.x + dx / len * 6, tip.y + dy / len * 6),
                                  ImVec2(tip.x - dy / len * 4, tip.y + dx / len * 4),
                                  ImVec2(tip.x + dy / len * 4, tip.y - dx / len * 4), kPlannerColor);
        }
    }

    // The inspected light ray.
    if (ctx.settings.rayInspector && ui.inspectorOpen && ui.inspectValid) {
        std::vector<ImVec2> pts;
        for (const auto& p : ui.inspect.path) pts.push_back(projectPoint(p));
        if (pts.size() > 1) dl->AddPolyline(pts.data(), static_cast<int>(pts.size()), kRayColor, 0, 1.5f);
    }

    if (nearest) {
        char label[96];
        std::snprintf(label, sizeof(label), "%s  r = %.2f M", nearest->name.c_str(), nearest->state.r);
        dl->AddText(ImVec2(nearestPos.x + 8, nearestPos.y - 8), IM_COL32(230, 230, 240, 255), label);
    }

    // Camera
    const Camera& cam = ctx.camera;
    auto cp = kerr::toCartesian(a, cam.r, cam.theta, cam.phi);
    ImVec2 c = project(static_cast<float>(cp[0]), static_cast<float>(cp[1]), static_cast<float>(cp[2]));
    dl->AddCircle(c, 5.0f, IM_COL32(120, 200, 255, 255), 12, 2.0f);
    dl->AddText(ImVec2(c.x + 7, c.y - 7), IM_COL32(120, 200, 255, 255), "camera");

    dl->PopClipRect();
    const char* caption = side ? "side view: spin axis up"
                        : ui.plannerOpen ? "top view: drag from a point to throw a body"
                                         : "top view: spin axis toward you";
    dl->AddText(ImVec2(origin.x + 6, maxCorner.y - ImGui::GetFontSize() - 4), IM_COL32(150, 150, 160, 255), caption);
    ImGui::End();
}

// --- Ray inspector -----------------------------------------------------------------

void inspectorWindow(UiContext& ctx) {
    UiState& ui = ctx.state;
    if (!ctx.settings.rayInspector) ui.inspectorOpen = false;  // switched off: forget the last result
    if (!ui.inspectorOpen || !ui.inspectValid) return;
    const bool insideHorizon = ctx.camera.insideHorizon(ctx.settings.spin);

    // Mark the inspected pixel.
    ImDrawList* fg = ImGui::GetForegroundDrawList();
    const ImVec2 p(ui.inspectX, ui.inspectY);
    fg->AddCircle(p, 9.0f, kRayColor, 20, 1.5f);
    fg->AddLine(ImVec2(p.x - 14, p.y), ImVec2(p.x - 5, p.y), kRayColor, 1.5f);
    fg->AddLine(ImVec2(p.x + 5, p.y), ImVec2(p.x + 14, p.y), kRayColor, 1.5f);
    fg->AddLine(ImVec2(p.x, p.y - 14), ImVec2(p.x, p.y - 5), kRayColor, 1.5f);
    fg->AddLine(ImVec2(p.x, p.y + 5), ImVec2(p.x, p.y + 14), kRayColor, 1.5f);

    const trajectory::RayReport& r = ui.inspect;
    ImGui::SetNextWindowSize(ImVec2(25 * ImGui::GetFontSize(), 0), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 30 * ImGui::GetFontSize(),
                                   ImGui::GetIO().DisplaySize.y - 16 * ImGui::GetFontSize()),
                            ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Ray inspector", &ui.inspectorOpen)) {
        ImGui::TextWrapped("The light reaching this pixel came from %s.", trajectory::rayEndName(r.end));
        switch (r.end) {
            case trajectory::RayEnd::Disk:
            case trajectory::RayEnd::Body: {
                const char* what = r.end == trajectory::RayEnd::Disk ? "disk gas" : ui.inspectBody.c_str();
                ImGui::TextWrapped("Emitted by %s at r = %.2f M. Frequency shift g = %.3f: it arrives %s (%s).",
                                   what, r.rEnd, r.g, r.g >= 1.0 ? "bluer and brighter" : "redder and dimmer",
                                   r.g >= 1.0 ? "moving toward you outweighs gravity"
                                              : "gravity and/or receding motion");
                break;
            }
            case trajectory::RayEnd::Horizon:
                if (insideHorizon) {
                    ImGui::TextWrapped("Traced backwards it never left the hole's interior: it fell in from the far "
                                       "side of the horizon, where nothing shines.");
                } else {
                    ImGui::TextWrapped("Traced backwards it ends at the horizon: no light can come from there. "
                                       "This pixel is part of the shadow.");
                }
                break;
            case trajectory::RayEnd::Sky: {
                double lon = std::atan2(r.skyDirection[1], r.skyDirection[0]) / kDeg;
                double lat = std::asin(std::clamp(r.skyDirection[2], -1.0, 1.0)) / kDeg;
                ImGui::TextWrapped("It left the sky from longitude %.1f, latitude %.1f degrees.", lon, lat);
                break;
            }
            case trajectory::RayEnd::Unfinished:
                break;
        }
        ImGui::Separator();
        ImGui::Text("Went around the hole %.2f times", r.turns);
        ImGui::Text("Closest approach  r = %.3f M", r.rMin);
        ImGui::Text("Travel time       %.1f M", r.delay);
        ImGui::Text("Impact parameter  b = %.3f M", r.b);
        ImGui::Text("Carter constant   q = %.3f M^2", r.q);
        ImGui::TextDisabled("Its path is drawn in yellow on the orbit map.");
        ImGui::TextDisabled("Ctrl+click or middle-click another pixel to inspect it.");
    }
    ImGui::End();
}

void helpWindow(UiState& ui) {
    if (!ui.showHelp) return;
    ImGui::SetNextWindowSize(ImVec2(460, 0), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Controls", &ui.showHelp)) {
        ImGui::TextUnformatted(
            "Left drag      orbit the camera around the hole\n"
            "Right drag     look around\n"
            "Scroll         move closer / farther\n"
            "W / S          move closer / farther\n"
            "A / D, Q / E   orbit sideways / up and down\n"
            "Ctrl+click     inspect the light reaching a pixel\n"
            "  (or middle-click)\n"
            "I              turn the ray inspector on or off\n"
            "Enter          launch the planned body\n"
            "Space          pause or resume orbits\n"
            "R              reset the camera\n"
            "G              toggle the lensing grid\n"
            "Tab            hide or show the interface\n"
            "F12            save a screenshot\n"
            "Esc            stop recording a video\n"
            "H              this window");
        ImGui::Separator();
        ImGui::TextWrapped(
            "Every pixel follows a light ray backwards through the curved spacetime of a rotating "
            "(Kerr) black hole. The dark region is the shadow: rays that end up in the horizon. "
            "The thin bright ring around it is light that circled the hole before escaping. The "
            "disk on the far side of the hole appears above and below it, bent over the top by "
            "gravitational lensing.");
    }
    ImGui::End();
}

}  // namespace

bool plannedLaunch(const UiState& ui, double spin, kerr::Constants& c, kerr::State& s) {
    double v[3];
    if (!plannedVelocity(ui, spin, v)) return false;
    kerr::launchFromZamo(spin, ui.launchR, ui.launchPolarDeg * kDeg, ui.launchAzimuthDeg * kDeg, v, 1.0, c, s);
    return true;
}

void drawUi(UiContext& ctx) {
    UiState& ui = ctx.state;
    helpWindow(ui);
    inspectorWindow(ctx);
    if (!ui.showPanel) {
        ui.plannerOpen = false;
        if (ui.showOrbitMap) orbitMap(ctx);
        return;
    }

    const float em = ImGui::GetFontSize();
    ImGui::SetNextWindowPos(ImVec2(12, 12), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(24.0f * em, ImGui::GetIO().DisplaySize.y - 24), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Black Hole Simulator")) {
        ImGui::PushItemWidth(-9.5f * em);  // leave room for labels
        statsSection(ctx);
        if (ImGui::CollapsingHeader("Black hole", ImGuiTreeNodeFlags_DefaultOpen)) blackHoleSection(ctx);
        if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) cameraSection(ctx);
        if (ImGui::CollapsingHeader("Accretion disk and jets", ImGuiTreeNodeFlags_DefaultOpen)) diskSection(ctx);
        ui.plannerOpen = false;
        if (ImGui::CollapsingHeader("Objects and orbits", ImGuiTreeNodeFlags_DefaultOpen)) objectsSection(ctx);
        if (ImGui::CollapsingHeader("Record video")) videoSection(ctx);
        if (ImGui::CollapsingHeader("Sky")) skySection(ctx);
        if (ImGui::CollapsingHeader("Rendering")) renderSection(ctx);
        ImGui::Separator();
        ImGui::Checkbox("Orbit map", &ui.showOrbitMap);
        ImGui::SameLine();
        ImGui::Checkbox("Controls", &ui.showHelp);
        ImGui::SameLine();
        ImGui::TextDisabled("Tab hides the interface");

        std::string errors = ctx.renderer.shaderErrors();
        if (!errors.empty()) {
            ImGui::SeparatorText("Shader errors");
            ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(255, 110, 110, 255));
            ImGui::TextWrapped("%s", errors.c_str());
            ImGui::PopStyleColor();
        }
        ImGui::PopItemWidth();
    }
    ImGui::End();
    if (ui.showOrbitMap) orbitMap(ctx);
}
