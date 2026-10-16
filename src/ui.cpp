#include "ui.h"

#include "blackbody.h"
#include "central.h"
#include "kerr.h"
#include "renderer.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
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

std::string printf1(const char* fmt, double v) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), fmt, v);
    return buf;
}

// A length or time for the readouts: in M, and in physical units when they are on.
std::string lengthText(const Settings& s, double lengthM) {
    std::string out = printf1("%.4g M", lengthM);
    if (s.physicalUnits) out += " (" + central::formatLength(lengthM, s.central.massSolar) + ")";
    return out;
}

std::string timeText(const Settings& s, double timeM) {
    std::string out = printf1("%.4g M", timeM);
    if (s.physicalUnits) out += " (" + central::formatTime(timeM, s.central.massSolar) + ")";
    return out;
}

double surfaceOf(const Settings& s) {
    return s.central.hasHorizon() ? 0.0 : s.central.radius;
}

// The custom launch's velocity in the ZAMO frame. False if it is not possible here.
bool plannedVelocity(const UiState& ui, double spin, double scale, double surface, double v[3]) {
    const double r = ui.launchR * scale;
    if (r < central::boundaryRadius(spin, surface) + 0.2 * scale) return false;
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

// A slider for a length kept in M but laid out in scene units: the range is given in
// scene units, the value shown in M.
bool sceneSlider(const char* label, float& sceneValue, float lo, float hi, double scale, bool logarithmic) {
    double value = sceneValue * scale;
    const double vMin = lo * scale, vMax = hi * scale;
    if (!ImGui::SliderScalar(label, ImGuiDataType_Double, &value, &vMin, &vMax, "%.4g M",
                             logarithmic ? ImGuiSliderFlags_Logarithmic : 0)) {
        return false;
    }
    sceneValue = static_cast<float>(value / scale);
    return true;
}

// --- Sections ------------------------------------------------------------------

void statsSection(UiContext& ctx) {
    const Renderer& r = ctx.renderer;
    ImGui::Text("%.0f fps   GPU %.1f ms   trace %dx%d", ctx.fps, r.gpuMilliseconds(), r.traceWidth(),
                r.traceHeight());
    ImGui::TextDisabled("%s", ctx.gpuName);
    ImGui::Text("t = %s", timeText(ctx.settings, ctx.bodies.time).c_str());
    ImGui::SameLine();
    ImGui::TextDisabled("(coordinate time of a distant observer)");
}

// The live numbers for the central object (surface physics, accretion, scales).
void centralReadouts(UiContext& ctx) {
    const Settings& s = ctx.settings;
    const central::Object& o = s.central;
    const double a = s.spin;
    if (!ImGui::BeginTable("readouts", 2, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) return;
    auto row = [](const char* name, const std::string& value, const char* tip) {
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::TextUnformatted(name);
        ImGui::SetItemTooltip("%s", tip);
        ImGui::TableNextColumn();
        ImGui::TextUnformatted(value.c_str());
    };
    if (!o.hasHorizon()) {
        const double R = o.radius;
        const double z = central::surfaceRedshift(R);
        const double gSI = central::accelerationSI(central::surfaceGravity(R), o.massSolar);
        row("Radius", printf1("%.5g M", R) + " = " + central::formatLength(R, o.massSolar),
            "In units of M = G M / c^2. The Sun is 471,000 M, a white dwarf ~10,000 M, a neutron star ~6 M.");
        row("Size", central::sizeComparison(central::lengthKm(R, o.massSolar)),
            "Its real size. Each object is framed to fill the view; the scale bar at the bottom of the view "
            "shows how big things are, and bodies around it are drawn at their real sizes.");
        row("Compactness M/R", printf1("%.3g", 1.0 / R) + printf1("   (%.3g x Buchdahl)", R / central::kBuchdahlRadius),
            "No static star can be smaller than the Buchdahl limit, R = 9/4 M: below it no pressure holds it up.");
        row("Surface redshift", printf1("1 + z = %.6g", 1.0 + z) + printf1("  (%.3g km/s)", z * central::kSpeedOfLightKmS),
            "Light climbing out of the surface loses energy: 1 + z = (1 - 2M/R)^(-1/2). Measured for the Sun "
            "(0.6 km/s) and for white dwarfs such as Sirius B (80 km/s).");
        row("Surface gravity", printf1("%.3g g", gSI / 9.80665) + printf1("  (%.3g m/s^2)", gSI),
            "Felt standing on the surface: M / (R^2 sqrt(1 - 2M/R)), which grows without limit near R = 2M.");
        row("Escape speed", printf1("%.4g c", central::escapeSpeed(R)), "sqrt(2M/R), measured at the surface.");
        row("Surface clocks", printf1("%.6g x a distant clock", central::surfaceClockRate(R)), "sqrt(1 - 2M/R).");
        static double cachedR = -1.0, cachedFraction = 0.0;
        if (R != cachedR) {
            cachedR = R;
            cachedFraction = central::visibleFraction(R);
        }
        row("Seen at once", printf1("%.1f%% of the surface", 100.0 * cachedFraction),
            "Light bending lets you see past the limb: below R = 3.52 M the whole surface is visible at once.");
        // Strong-gravity features only matter for compact objects: for a star or planet
        // they lie thousands of times deeper than its surface.
        if (R < 30.0) {
            row("Photon ring", R < 3.0 ? "outside the surface" : "inside the star (none)",
                "Light can circle at r = 3M. Only a star smaller than that shows the photon ring.");
            const double isco = kerr::iscoRadius(a, true);
            row("ISCO", R < isco ? lengthText(s, isco) + " (outside: disk stops there)" : std::string("inside the star"),
                "Innermost stable circular orbit. A neutron star (R ~ 5.8 M) is smaller than its ISCO.");
        }
        const central::SurfaceAccretion acc = central::surfaceEfficiency(a, R, o.surfaceOmega());
        row("Accretion: disk", printf1("%.4g%% of rest mass", 100.0 * acc.disk),
            "Energy released by gas spiralling through the disk to its inner edge, max(ISCO, R).");
        row("Accretion: surface", printf1("%.4g%%", 100.0 * acc.boundaryLayer),
            "Released when the gas settles onto the surface (the boundary layer). A black hole swallows this.");
        row("Accretion: total", printf1("%.4g%%", 100.0 * acc.total),
            "1 - E_surface. For a neutron star more than a black hole of the same spin, because nothing falls "
            "through a horizon unseen.");
        const bool magnetic = o.kind == central::Kind::NeutronStar || o.kind == central::Kind::WhiteDwarf;
        if (magnetic && o.magneticField > 0.0) {
            const double rA = central::alfvenRadiusKm(o.magneticField, central::lengthKm(R, o.massSolar), o.massSolar,
                                                      s.accretionRate * central::eddingtonRate(o.massSolar)) /
                              (central::kKmPerSolarMass * o.massSolar);
            row("Alfven radius", lengthText(s, rA) + (0.5 * rA > R ? "" : " (below the surface)"),
                "Where the magnetic field takes over from the disk: r_A = (mu^4 / (2 G M Mdot^2))^(1/7). The disk "
                "stops near half of it; inside, gas flows along field lines onto the magnetic poles.");
        }
    } else {
        row("Size", "horizon " + central::sizeComparison(central::lengthKm(kerr::horizonRadius(a), o.massSolar)),
            "The event horizon's real size for this mass: 3 km across per solar mass without spin.");
        row("Accretion efficiency", printf1("%.4g%% of rest mass", 100.0 * central::blackHoleEfficiency(a)),
            "1 - E_isco: the energy a thin disk releases before the gas crosses the horizon. 5.7% without spin, "
            "32% at a = 0.998.");
        const double rStar = central::kSunRadiusKm / (central::kKmPerSolarMass * o.massSolar);
        const double rt = central::tidalRadius(rStar, 1.0 / o.massSolar);
        const double rmb = kerr::marginallyBoundRadius(a, true);
        row("Sun-like star torn at", rt < rmb ? std::string("(swallowed whole)") : lengthText(s, rt),
            "Its tidal radius, r_t = R* (M / m*)^(1/3). Holes heavier than the Hills mass (1.1e8 Suns for a Sun-like "
            "star) swallow such stars before the tides can tear them.");
    }
    row("1 M is", central::formatLength(1.0, o.massSolar) + ", " + central::formatTime(1.0, o.massSolar),
        "The units of length and time, G M / c^2 and G M / c^3, for this mass.");
    ImGui::EndTable();
}

void radiiTable(UiContext& ctx) {
    const Settings& s = ctx.settings;
    const double a = s.spin;
    const double surface = surfaceOf(s);
    if (!ImGui::BeginTable("radii", 3, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) return;
    ImGui::TableSetupColumn("Radius (M)");
    ImGui::TableSetupColumn("Prograde");
    ImGui::TableSetupColumn("Retrograde");
    ImGui::TableHeadersRow();
    auto value = [&](double r) {
        if (r < 0.0) return;
        if (surface > 0.0 && r < surface) ImGui::TextDisabled("%.3f (inside)", r);
        else ImGui::Text("%.3f", r);
    };
    auto row = [&](const char* name, double pro, double retro, const char* tip) {
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::TextUnformatted(name);
        ImGui::SetItemTooltip("%s", tip);
        ImGui::TableNextColumn();
        value(pro);
        ImGui::TableNextColumn();
        value(retro);
    };
    if (surface > 0.0) {
        row("Surface", surface, -1, "Rays, bodies, gas and the camera stop here.");
    } else {
        row("Event horizon", kerr::horizonRadius(a), -1, "Nothing, not even light, returns from inside.");
        row("Ergosphere (equator)", kerr::ergosphereRadius(a, 0.5 * kerr::kPi), -1,
            "Inside, frame dragging is so strong that nothing can stay still relative to the stars.");
    }
    row("Photon orbit", kerr::photonOrbitRadius(a, true), kerr::photonOrbitRadius(a, false),
        "Light can circle the hole here (unstably). It forms the bright photon ring.");
    row("ISCO", kerr::iscoRadius(a, true), kerr::iscoRadius(a, false),
        "Innermost stable circular orbit: where the accretion disk ends.");
    ImGui::EndTable();
}

void centralSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    UiState& ui = ctx.state;
    central::Object& o = s.central;
    const auto& list = central::catalog();

    // The catalogue, one click each.
    const float buttonWidth = (ImGui::GetContentRegionAvail().x - 2.0f * ImGui::GetStyle().ItemSpacing.x) / 3.0f;
    for (size_t k = 0; k < list.size(); ++k) {
        if (k % 3 != 0) ImGui::SameLine();
        if (ImGui::Button(list[k].label, ImVec2(buttonWidth, 0))) ctx.centralRequested = static_cast<int>(k);
        ImGui::SetItemTooltip("%s", list[k].description);
    }
    ImGui::TextDisabled("Famous black holes:");
    for (int i = 1; i < SceneCount; ++i) {
        if (i > 1) ImGui::SameLine();
        if (ImGui::SmallButton(kSceneNames[i])) ctx.sceneRequested = i;
        ImGui::SetItemTooltip("%s", kSceneDescriptions[i]);
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Reset")) ctx.sceneRequested = SceneDefault;

    ImGui::SeparatorText((o.name + "  (" + central::kindName(o.kind) + ")").c_str());
    const bool neutron = o.kind == central::Kind::NeutronStar, dwarf = o.kind == central::Kind::WhiteDwarf;

    // Mass: for a body, the radius in km stays (or follows the mass-radius relation).
    {
        double mass = o.massSolar;
        const double lo = 1e-7, hi = 3e10;
        if (ImGui::SliderScalar("Mass", ImGuiDataType_Double, &mass, &lo, &hi, "%.4g Msun", ImGuiSliderFlags_Logarithmic)) {
            const double km = central::lengthKm(o.radius, o.massSolar);
            o.massSolar = mass;
            if (!o.hasHorizon()) {
                double followed = 0.0;
                if (s.radiusFollowsMass && dwarf) followed = central::whiteDwarfRadiusKm(mass);
                if (s.radiusFollowsMass && neutron) followed = central::neutronStarRadiusKm(mass);
                o.radius = (followed > 0.0 ? followed : km) / (central::kKmPerSolarMass * mass);
            }
            ui.compactnessOrigin = -1;
        }
        help("For a black hole the mass only sets the units: everything scales with it. A star or planet keeps "
             "its size in km (or follows the mass-radius relation of white dwarfs and neutron stars), so its R/M "
             "changes. White dwarfs cannot exceed 1.44 solar masses, neutron stars ~1.64 (this equation of state).");
    }
    if (!o.hasHorizon()) {
        if (ui.radiusInKm) {
            double km = central::lengthKm(o.radius, o.massSolar);
            const double lo = 1.0, hi = 1e9;
            if (ImGui::SliderScalar("Radius", ImGuiDataType_Double, &km, &lo, &hi, "%.4g km", ImGuiSliderFlags_Logarithmic)) {
                o.radius = km / (central::kKmPerSolarMass * o.massSolar);
                ui.compactnessOrigin = -1;
            }
        } else {
            double R = o.radius;
            const double lo = 2.0, hi = 2e9;
            if (ImGui::SliderScalar("Radius", ImGuiDataType_Double, &R, &lo, &hi, "%.5g M", ImGuiSliderFlags_Logarithmic)) {
                o.radius = R;
                ui.compactnessOrigin = -1;
            }
        }
        ImGui::Checkbox("Radius in km", &ui.radiusInKm);
        if (neutron || dwarf) {
            ImGui::SameLine();
            ImGui::Checkbox("Radius follows mass", &s.radiusFollowsMass);
            help("White dwarfs: Nauenberg's relation, R = 0.0112 R_sun sqrt((M/M_Ch)^(-2/3) - (M/M_Ch)^(2/3)): "
                 "heavier ones are smaller. Neutron stars: the Tolman-Oppenheimer-Volkoff equations for a "
                 "Gamma = 2 polytrope (14 km at 1.4 solar masses).");
        }
    }

    // Spin.
    {
        double chi = o.chi;
        const double lo = 0.0, hi = o.hasHorizon() ? 0.999 : o.maxChi();
        if (ImGui::SliderScalar(o.hasHorizon() ? "Spin a/M" : "Spin J/M^2", ImGuiDataType_Double, &chi, &lo, &hi, "%.4g",
                                hi > 5.0 ? ImGuiSliderFlags_Logarithmic : 0)) {
            o.chi = chi;
        }
        if (o.hasHorizon()) {
            help("Angular momentum per unit mass. 0 is a non-rotating (Schwarzschild) hole; values near 1 are extreme "
                 "Kerr holes whose rotation drags space around them. Most astrophysical black holes are thought to "
                 "spin quickly.");
        } else {
            help("c J / (G M^2). Above 1 is normal for planets and ordinary stars (the Earth's is ~890); only black "
                 "holes cannot exceed 1. The slider stops at the mass-shedding limit, where the equator would fly "
                 "apart (Omega ~ 0.65 sqrt(M / R^3)). Outside, the Kerr metric with a = min(J/M^2, 0.99) is used: "
                 "right in mass and angular momentum, approximate in the quadrupole.");
            const std::string period =
                central::formatTime(2.0 * kerr::kPi / std::max(o.surfaceOmega(), 1e-300), o.massSolar);
            if (o.radius < 30.0) {
                ImGui::TextDisabled("Turns once every %s; Kerr spin a = %.3f", period.c_str(), s.spin);
            } else {
                // Frame dragging at the surface, omega ~ 2J/R^3, against the surface's own spin.
                ImGui::TextDisabled("Turns once every %s; frame dragging %.2g of that", period.c_str(),
                                    2.0 * o.chi / (o.radius * o.radius * o.radius) / std::max(o.surfaceOmega(), 1e-300));
            }
            if (o.chi > o.maxChi() * 1.0001) {
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.65f, 0.3f, 1.0f));
                ImGui::TextWrapped("Squeezed, it would spin faster than its surface can (mass shedding): the surface "
                                   "turns at the limit, and J is kept for a collapse.");
                ImGui::PopStyleColor();
            }
        }
    }

    if (!o.hasHorizon()) {
        float T = static_cast<float>(o.surfaceTemperature);
        if (ImGui::SliderFloat("Surface temperature", &T, 100.0f, 1e8f, "%.0f K", ImGuiSliderFlags_Logarithmic)) {
            o.surfaceTemperature = T;
        }
        if (T > 30000.0f) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
            ImGui::TextWrapped("Drawn as %.0f K: hotter surfaces shine mostly in UV and X-rays, so their color "
                               "temperature is compressed.", central::displayTemperature(T));
            ImGui::PopStyleColor();
        }
        if (neutron || dwarf) {
            double B = std::max(o.magneticField, 1.0);
            const double lo = 1.0, hi = 1e16;
            if (ImGui::SliderScalar("Magnetic field", ImGuiDataType_Double, &B, &lo, &hi, "%.3g G",
                                    ImGuiSliderFlags_Logarithmic)) {
                o.magneticField = B;
            }
            ImGui::SliderFloat("Accretion rate", &s.accretionRate, 1e-4f, 1.0f, "%.3g Eddington", ImGuiSliderFlags_Logarithmic);
            help("How fast the disk feeds the star, in units of the Eddington rate (where radiation pressure would "
                 "balance gravity). With the field, it sets where the magnetosphere stops the disk.");
            ImGui::Checkbox("Magnetosphere", &s.magnetosphere);
            help("A magnetized star stops the disk near half the Alfven radius; inside, the gas follows the dipole "
                 "field lines r = L sin^2(theta) down onto the magnetic poles (accretion columns), heating them.");
            ImGui::SameLine();
            if (neutron) ImGui::Checkbox("Pulsar beams", &o.pulsar);
            ImGui::SameLine();
            if (ImGui::SmallButton(neutron ? "X-ray burst" : "Nova")) ctx.burstRequested = true;
            ImGui::SetItemTooltip(neutron ? "Helium accreted on a neutron star ignites and burns across the whole surface "
                                            "in seconds (a Type I X-ray burst)."
                                          : "Hydrogen accreted on a white dwarf ignites in a thermonuclear runaway.");
            if (neutron) {
                ImGui::Checkbox("Field lines", &s.fieldLines);
                help("The dipole's closed field lines, r = L sin^2(theta) about the magnetic axis, out to the light "
                     "cylinder c / Omega or the disk's edge, turning with the star and bent by its gravity like "
                     "everything else. A magnetar's are twisted (~1 radian) and wrapped in a glow of electron-"
                     "positron pairs. False color: this plasma shines in X-rays and radio, not visible light.");
                ImGui::SameLine();
                ImGui::Checkbox("Beams", &s.pulsarBeams);
                help("A pulsar's twin beams, from the open field lines over its magnetic poles, at the opening angle "
                     "the last open field line gives (28 degrees for the Crab, ~4 for a slow magnetar). Drawn as if "
                     "they lit a faint haze; light far out left when the star had turned less, so fast pulsars' "
                     "beams curve into a spiral. A real beam is seen only as a flash as it sweeps over you.");
                ImGui::SliderFloat("Magnetosphere glow", &s.magnetosphereGlow, 0.05f, 10.0f, "%.2f",
                                   ImGuiSliderFlags_Logarithmic);
            }
        }
        ImGui::Checkbox("Surface detail", &s.surfaceDetail);
        help("Granulation (convection cells sized from the pressure scale height, boiling at their real "
             "lifetime), starspots in their real latitudes and groups, faculae, flares on red dwarfs, hot "
             "polar caps and flares on neutron stars, impact hot spots. Fine detail fades in as you get close.");
        ImGui::SameLine();
        ImGui::Checkbox("Corona / atmosphere", &s.corona);
        help("A star's corona (and pink chromosphere), or a planet's air. The air is drawn at true scale: sunlight "
             "scattered by molecules (blue), aerosols and ozone, 100 km thick on an Earth-sized planet. The "
             "corona has the real Sun's profile and streamers, brightened (below) since the real one is a "
             "millionth as bright as the disk.");
        if (o.kind == central::Kind::RockyPlanet || o.kind == central::Kind::GasGiant) {
            ImGui::Checkbox("Spacecraft maps", &s.planetMaps);
            help("The Earth from NASA's Blue Marble (with the sea floor), its cloud cover and the Black Marble's city "
                 "lights; Jupiter from Cassini's global map. The light, the air, the glint and the night side are "
                 "still computed. Untick for the procedural planet (used for any planet that is not these two).");
        }
        if (o.kind == central::Kind::Star || dwarf || neutron) {
            ImGui::Checkbox("Film colors", &s.filmColor);
            help("A blackbody's true color is pale: the Sun is white, a 3,000 K red dwarf peach. Films and colorized "
                 "solar images show stars warmer and deeper: the Sun golden with an orange limb, red dwarfs red. "
                 "Only the hue changes; the brightness (limb darkening, spots, granulation, flares) is the physics. "
                 "Untick for the true color.");
            if (dwarf && s.diskEnabled) {
                ImGui::SameLine();
                ImGui::Checkbox("Companion", &s.companionStar);
                help("The red dwarf feeding the disk (0.15 solar masses), filling its Roche lobe 57 radii away, and "
                     "the stream of its gas falling from the L1 point along its ballistic path to the disk's rim.");
            }
        }
        if (o.kind == central::Kind::Star && s.corona) {
            ImGui::SliderFloat("Corona brightness", &s.coronaBoost, 1.0f, 1e6f, "x%.0f real", ImGuiSliderFlags_Logarithmic);
            help("x1 is the real Sun's corona, a millionth of the disk's brightness, seen only when the disk is "
                 "hidden (a total eclipse).");
        }
        ImGui::SliderFloat("Surface brightness", &s.surfaceBrightness, 0.05f, 20.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
    }

    // Compactness: keep the mass, change the radius.
    if (!o.hasHorizon()) {
        double C = 1.0 / o.radius;
        const double cMin = std::min(2e-6, 1.0 / o.radius), cMax = 0.5;
        if (ImGui::SliderScalar("Compactness M/R", ImGuiDataType_Double, &C, &cMin, &cMax, "%.3g",
                                ImGuiSliderFlags_Logarithmic)) {
            if (ui.compactnessOrigin < 0) {
                ui.compactnessOrigin = static_cast<int>(o.kind);
                ui.compactnessTemperature = o.surfaceTemperature;
            }
            o.radius = 1.0 / C;
            if (s.autoStage) {
                const central::Kind stage =
                    central::stageForCompactness(o.radius, static_cast<central::Kind>(ui.compactnessOrigin));
                if (central::becomeStage(o, stage, ui.compactnessTemperature)) {
                    if (stage == central::Kind::WhiteDwarf) {
                        ui.centralStatus = "As dense as a white dwarf: a teaspoon weighs tonnes, and degenerate electrons "
                                           "hold it up (real white dwarfs have R ~ 1,000-30,000 M).";
                    } else if (stage == central::Kind::NeutronStar) {
                        ui.centralStatus = "Neutron-star density: electrons are squeezed into the protons, and degenerate "
                                           "neutrons hold it up. (No stable star exists between ~10 and ~1,000 M: there "
                                           "a shrinking star collapses.)";
                    } else {
                        ui.centralStatus = "Swollen back to an ordinary star.";
                    }
                }
            }
        }
        help("Keeps the mass and shrinks (or swells) the radius. Outside the surface nothing changes (Birkhoff's "
             "theorem): orbits and light bending stay those of a black hole of the same mass. What changes is how "
             "deep you can go: below R = 6M the disk ends at the ISCO, below 3.5M the whole surface shows, below 3M "
             "the photon ring appears, and below the Buchdahl limit 9/4 M nothing can stand: it collapses into a "
             "black hole with the same mass and spin.");
        ImGui::Checkbox("Becomes a white dwarf, then a neutron star", &s.autoStage);
    }
    ImGui::Checkbox("Physical units", &s.physicalUnits);
    help("Readouts also in km, AU, seconds and hours, for the mass above. With a black hole, the Tidal disruption "
         "preset then uses a Sun-like star at its true size.");

    if (!ui.centralStatus.empty()) ImGui::TextWrapped("%s", ui.centralStatus.c_str());
    if (!o.hasHorizon()) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.6f, 0.85f, 1.0f, 1.0f));
        ImGui::TextWrapped(o.radius < 30.0
                               ? "Gravity outside depends only on the mass and spin (Birkhoff's theorem): out here a "
                                 "black hole of the same mass would pull the same way. What differs is the surface, "
                                 "and how close you can get."
                               : "Its gravity is ordinary Newtonian gravity to a few parts in a million (M/R): light "
                                 "bends by microradians, clocks barely slow. Relativity only shows if you squeeze it "
                                 "(Compactness, above).");
        ImGui::PopStyleColor();
    }
    centralReadouts(ctx);
    if (o.hasHorizon() || o.radius < 30.0) radiiTable(ctx);
}

void cameraSection(UiContext& ctx) {
    Camera& cam = ctx.camera;
    Settings& s = ctx.settings;
    UiState& ui = ctx.state;
    const double S = ctx.scale;
    const double minR = cam.minDistance(s.spin);

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
    ImGui::SetItemTooltip("Drop from rest straight toward the center. Toward a black hole: watch the shadow swell, "
                          "then cross the horizon. Toward a star or planet: you land on its surface.");

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
            ImGui::Text("Speed %.3f c   gamma %.2f   r = %s", speed, 1.0f / std::sqrt(std::max(1.0f - speed * speed, 1e-6f)),
                        lengthText(s, cam.r).c_str());
            ImGui::Text("Your clock: %.4g M (runs at %.0f%% of a distant clock)", cam.properTime,
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
        ImGui::TextDisabled("Clocks here run at %.4g%% of the rate far away.", 100.0 * cam.clockRate(s.spin));
    }
    if (!ui.cameraStatus.empty()) ImGui::TextWrapped("%s", ui.cameraStatus.c_str());
    if (motion == CameraMotion::Hover) {
        double r = cam.r;
        const double rMax = 500.0 * S;
        if (ImGui::SliderScalar("Distance", ImGuiDataType_Double, &r, &minR, &rMax, "%.4g M",
                                ImGuiSliderFlags_Logarithmic)) {
            cam.r = r;
        }
        if (s.physicalUnits) ImGui::TextDisabled("= %s", central::formatLength(cam.r, s.central.massSolar).c_str());
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
                cam.r = std::max(rr * S, minR);
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
    const double S = ctx.scale;
    ImGui::Checkbox("Show disk", &s.diskEnabled);
    ImGui::SameLine();
    if (ImGui::RadioButton("Prograde", s.diskPrograde)) s.diskPrograde = true;
    ImGui::SameLine();
    if (ImGui::RadioButton("Retrograde", !s.diskPrograde)) s.diskPrograde = false;

    ImGui::Checkbox(s.central.hasHorizon() ? "Inner edge at ISCO" : "Inner edge at ISCO or surface", &s.diskInnerAtIsco);
    if (!s.diskInnerAtIsco) sceneSlider("Inner radius", s.diskInner, 1.0f, 15.0f, S, false);
    const bool hole = s.central.hasHorizon();
    sceneSlider("Outer radius", s.diskOuter, 4.0f, hole ? 60.0f : 200.0f, S, !hole);
    if (s.diskPhysicalTemperature > 0.0) {
        // Around a neutron star or white dwarf the temperature follows from the physics.
        ImGui::Text("Peak temperature %s K", printf1("%.3g", s.diskPhysicalTemperature).c_str());
        help("A thin disk's hottest ring, sigma T^4 = 3 G M Mdot / (8 pi r^3) (1 - sqrt(r_in / r)) at r = 49/36 "
             "r_in, set by the mass, the accretion rate and the inner edge: ~1e7 K (X-rays) around a neutron "
             "star, a few 1e4 K (ultraviolet) around a white dwarf. Drawn at a compressed color temperature, "
             "as hot surfaces are.");
        ImGui::SameLine();
        ImGui::TextDisabled("(drawn as %.0f K)", s.diskTemperature);
    } else {
        ImGui::SliderFloat("Peak temperature", &s.diskTemperature, 1500.0f, 40000.0f, "%.0f K",
                           ImGuiSliderFlags_Logarithmic);
        help("Sets the disk's color through the blackbody spectrum. Real disks around stellar black holes "
             "reach millions of kelvin; here it is kept in the visible range so the colors mean something.");
    }
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
    // Only what could really be found around this center.
    const central::Host host = central::hostOf(ctx.settings.central);
    for (const PresetInfo& p : presetList()) {
        if (p.object != objects || !presetAvailable(p, host)) continue;
        if (column++ % 3 != 0) ImGui::SameLine();
        if (ImGui::Button(p.label, ImVec2(buttonWidth, 0))) ctx.bodies.addPreset(ctx.settings.spin, p.preset);
        ImGui::SetItemTooltip("%s", p.description);
    }
    if (objects && column == 0) {
        ImGui::TextDisabled("Nothing is known to orbit this kind of object; try the orbits below.");
    }
    if (objects && host == central::Host::BlackHole) {
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
    const double S = ctx.scale, tS = S * std::sqrt(S);
    const double surface = surfaceOf(ctx.settings);
    char key[320];
    std::snprintf(key, sizeof(key), "%.6f %.4f %.3f %.3f %.4f %.3f %.4f %d %.6g %.6g", a, ui.launchR, ui.launchPolarDeg,
                  ui.launchAzimuthDeg, ui.launchSpeed, ui.launchTiltDeg, ui.launchRadial, ui.launchPrograde ? 1 : 0, S,
                  surface);
    if (ui.plannerKey == key) return;
    ui.plannerKey = key;
    ui.plannerValid = plannedLaunch(ui, a, S, surface, ui.plannerConstants, ui.plannerState);
    if (!ui.plannerValid) return;
    ui.plannerSpin = a;
    ui.plannerAnalysis = trajectory::classify(a, ui.plannerConstants, ui.plannerState, surface);
    ui.plannerPath = trajectory::predict(a, ui.plannerConstants, ui.plannerState, 3000.0 * tS, 0.25 * S, 8000, surface,
                                         2000.0 * S);
}

// The effective potential V(r): the body can only be where its energy E is above it.
void potentialPlot(const UiState& ui, double a, double S, double surface) {
    const kerr::Constants& c = ui.plannerConstants;
    const trajectory::Analysis& an = ui.plannerAnalysis;
    const double rh = central::boundaryRadius(a, surface);
    const double r0 = ui.plannerState.r;
    double rMax = std::max(40.0 * S, 1.6 * r0);
    if (an.hasOuter) rMax = std::max(rMax, 1.3 * an.rOuter);
    rMax = std::min(rMax, 2000.0 * S);
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
    // Far from a big object the well is shallow: zoom in on it.
    const double span = std::max(std::min(vPeak, E + 0.25) - vMin, 1e-12);
    double yHi = std::max(E, std::min(vPeak, E + 0.25)) + 0.12 * span;
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
    // ISCO marker.
    const double isco = kerr::iscoRadius(a, true);
    if (isco < rMax && isco > rMin) {
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
    char label[96];
    std::snprintf(label, sizeof(label), "%s ... r (log) ... %.4g M", surface > 0.0 ? "surface" : "horizon", rMax);
    const float labelWidth = ImGui::CalcTextSize(label).x;
    dl->AddText(ImVec2(corner.x - labelWidth - 4, origin.y + 2), IM_COL32(150, 150, 160, 255), label);
}

void plannerSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    BodySystem& sys = ctx.bodies;
    UiState& ui = ctx.state;
    const double a = s.spin;
    const double S = ctx.scale;
    const float fS = static_cast<float>(S);
    const double surface = surfaceOf(s);
    ui.plannerOpen = true;

    ImGui::TextWrapped("Set a start and a velocity, or drag on the orbit map (top view) from the start "
                       "point in the direction to throw. The verdict and the dashed path update live.");
    const int previousKind = ui.launchKind;
    ImGui::Combo("Type", &ui.launchKind,
                 "Star\0Rocky planet\0Gas giant\0Pulsar\0Star (tidally disrupted)\0White dwarf\0Neutron star\0"
                 "Small black hole\0Moon or rock\0");
    // Around a star or planet, bodies have real sizes: a typical one of the chosen kind.
    const bool realSizes = !s.central.hasHorizon();
    if (realSizes && (ui.launchKind != previousKind || ui.launchRadiusKm <= 0.0f)) {
        static const float typicalKm[] = {695700.0f, 6371.0f, 69911.0f, 12.0f, 695700.0f, 8700.0f, 12.0f, 0.0f, 1737.0f};
        ui.launchRadiusKm = typicalKm[std::clamp(ui.launchKind, 0, 8)];
    }
    sceneSlider("Start radius", ui.launchR, 2.0f, 80.0f, S, true);
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
    const bool hole = kind == BodyKind::MiniBlackHole;
    const bool rock = kind == BodyKind::Moon;
    if (!planet && !hole && !rock) {
        ImGui::SliderFloat("Temperature", &ui.launchTemperature, 2000.0f, 40000.0f, "%.0f K", ImGuiSliderFlags_Logarithmic);
    }
    if (!hole && realSizes) {
        ImGui::SliderFloat("Radius", &ui.launchRadiusKm, 0.01f, 1e7f, "%.4g km", ImGuiSliderFlags_Logarithmic);
        help("Its real radius. The Sun is 695,700 km, Jupiter 69,911, the Earth 6,371, the Moon 1,737, a white "
             "dwarf ~8,700, a neutron star ~12. Smaller than a pixel, it shows as a point of light.");
    } else if (!hole) {
        sceneSlider("Size", ui.launchRadius, 0.1f, 2.0f, S, false);
    }
    if (kind != BodyKind::TidalStar) {
        ImGui::SliderFloat("Mass", &ui.launchMass, 1e-7f, 0.05f, "%.2g M", ImGuiSliderFlags_Logarithmic);
        help("In units of the central mass. A Sun-like star around a million-Sun hole is 1e-6: its pull "
             "on anything else is negligible. Raise it (to ~0.01) and, with gravity between bodies on, it "
             "can hold moons and capture or fling other bodies. A small black hole's mass sets how strongly it "
             "lenses (4m/b) and the size of its shadow (3 sqrt(3) m).");
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
            case trajectory::Fate::HitsSurface: {
                const bool surfaceHit = an.fate == trajectory::Fate::HitsSurface;
                const char* into = surfaceHit ? "Hits the surface" : "Falls into the black hole";
                color = ImVec4(1.0f, 0.4f, 0.35f, 1.0f);
                if (path.fate == an.fate) {
                    std::snprintf(verdict, sizeof(verdict), "%s after %.4g M.", into, path.duration);
                } else if (an.hasOuter) {
                    std::snprintf(verdict, sizeof(verdict), "Swings out to r = %.4g M, then %s.", an.rOuter,
                                  surfaceHit ? "hits the surface" : "falls in");
                } else {
                    std::snprintf(verdict, sizeof(verdict), "%s.", into);
                }
                break;
            }
            case trajectory::Fate::Escapes:
                color = ImVec4(1.0f, 0.85f, 0.3f, 1.0f);
                std::snprintf(verdict, sizeof(verdict), "Escapes to infinity%s (E = %.4f).",
                              an.hasInner ? " after swinging past the center" : "", ui.plannerConstants.E);
                break;
            case trajectory::Fate::Circular:
                if (!an.stableCircular) color = ImVec4(1.0f, 0.65f, 0.3f, 1.0f);
                std::snprintf(verdict, sizeof(verdict), "%s", an.stableCircular
                                  ? "Stable circular orbit."
                                  : "Circular, but unstable: the slightest nudge sends it in or out.");
                break;
            case trajectory::Fate::Bound:
                std::snprintf(verdict, sizeof(verdict), "Bound orbit between r = %.4g and %.4g M.", an.rInner,
                              an.rOuter);
                break;
        }
        ImGui::PushStyleColor(ImGuiCol_Text, color);
        ImGui::TextWrapped("%s", verdict);
        ImGui::PopStyleColor();
        ImGui::TextDisabled("E = %.6f   L = %+.4g   Q = %.4g   tilt %.1f deg", ui.plannerConstants.E,
                            ui.plannerConstants.L, ui.plannerConstants.Q, an.inclinationDeg);
        if (an.fate == trajectory::Fate::Bound && path.periapses >= 2) {
            ImGui::TextDisabled("The closest point advances %.3g deg every orbit.", path.precessionDeg);
        }
        if (ui.plannerJointValid) {
            const trajectory::Path& joint = ui.plannerJointPath;
            char among[256];
            if (joint.fate == trajectory::Fate::Plunges || joint.fate == trajectory::Fate::HitsSurface) {
                std::snprintf(among, sizeof(among), "With the other bodies' pull: %s after %.4g M.",
                              joint.fate == trajectory::Fate::Plunges ? "falls in" : "hits the surface", joint.duration);
            } else if (joint.fate == trajectory::Fate::Escapes) {
                std::snprintf(among, sizeof(among), "With the other bodies' pull: escapes (gone after %.4g M).",
                              joint.duration);
            } else {
                std::snprintf(among, sizeof(among),
                              "With the other bodies' pull: still in orbit after %.4g M, between r = %.4g and %.4g M.",
                              joint.duration, joint.rMin, joint.rMax);
            }
            ImGui::TextWrapped("%s", among);
            help("The verdict above is exact for a body that only feels the central mass. This line comes from "
                 "simulating it together with the other bodies (gravity between bodies is on), which can "
                 "change its fate. The dashed path shows this simulation.");
        }
        potentialPlot(ui, a, S, surface);
        help("Effective potential: the least energy a body with this angular momentum needs to be at each "
             "radius. It can only go where its energy line (blue) is above the curve; the dots are the turning "
             "points. A hump the line clears means nothing stops it falling in.");
    } else {
        ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.35f, 1.0f), surface > 0.0 ? "Cannot start there: inside the surface."
                                                                          : "Cannot start there: too close to the horizon.");
    }

    auto finish = [&]() {
        Body& b = sys.bodies.back();
        b.temperature = planet || hole || rock ? 0.0f : ui.launchTemperature;
        b.radius = realSizes ? static_cast<float>(central::kmToM(ui.launchRadiusKm, s.central.massSolar))
                             : ui.launchRadius * fS;
        const float spinScale = 1.0f / (fS * std::sqrt(fS));
        if (kind == BodyKind::RockyPlanet) b.spinRate = 0.06f * spinScale;
        if (kind == BodyKind::GasGiant) b.spinRate = 0.1f * spinScale;
        if (kind == BodyKind::Pulsar) b.spinRate = 2.0f * spinScale;
        b.mass = ui.launchMass;
        if (hole) b.radius = static_cast<float>(3.0 * std::sqrt(3.0) * ui.launchMass);
        if (kind == BodyKind::TidalStar) {
            b.tidalRadius = static_cast<float>(std::max(15.0, 0.6 * ui.launchR) * S);
            b.mass = std::pow(b.radius / b.tidalRadius, 3.0);  // r_t = R (M / m)^(1/3)
        }
    };
    const bool enter = !ImGui::GetIO().WantTextInput && ImGui::IsKeyPressed(ImGuiKey_Enter, false);
    if (ImGui::Button("Launch (Enter)") || enter) {
        double v[3];
        if (plannedVelocity(ui, a, S, surface, v) &&
            sys.addFromZamo(a, ui.launchR * S, ui.launchPolarDeg * kDeg, ui.launchAzimuthDeg * kDeg, v, kind)) {
            finish();
        }
    }
    if (kind == BodyKind::TidalStar) {
        ImGui::SameLine();
        ImGui::TextDisabled("Tidal radius %.4g M: torn apart well inside it.", std::max(15.0, 0.6 * ui.launchR) * S);
    }
    ImGui::SeparatorText("Or by turning points (equatorial)");
    sceneSlider("Periapsis", ui.launchRp, 1.0f, 40.0f, S, false);
    sceneSlider("Apoapsis", ui.launchRa, 2.0f, 80.0f, S, false);
    if (ImGui::Button("Launch bound orbit")) {
        if (sys.addEquatorial(a, ui.launchRp * S, ui.launchRa * S, ui.launchPrograde, ui.launchAzimuthDeg * kDeg, kind)) {
            finish();
        }
    }
}

void objectsSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    BodySystem& sys = ctx.bodies;
    const double S = ctx.scale, tS = S * std::sqrt(S);

    ImGui::Checkbox("Paused", &s.paused);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-1);
    char speed[96];
    std::snprintf(speed, sizeof(speed), "%%.4g x   (%.3g M per second)", s.timeScale * tS);
    // Moons and planets take days to years to go round: around a star or planet time can
    // run far faster.
    const bool hole = s.central.hasHorizon();
    ImGui::SliderFloat("##timescale", &s.timeScale, 0.1f, hole ? 500.0f : 1e5f, speed, ImGuiSliderFlags_Logarithmic);
    if (s.physicalUnits || !hole) {
        ImGui::TextDisabled("Each second shows %s.", central::formatTime(s.timeScale * tS, s.central.massSolar).c_str());
    }
    ImGui::Checkbox("Light travel time", &s.lightDelay);
    help("Light takes time to reach you, so everything is seen where it was when its light set out: "
         "fast stars lag behind and stretch, lensed images of the same star show it at different "
         "moments, and a pulsar's lensed images flash out of step.");
    ImGui::Checkbox("Gravity between bodies", &s.bodyGravity);
    help("Bodies pull on each other (and on debris gas), each feeling the others' Newtonian pull in its "
         "own rest frame on top of the central mass's exact gravity. Most bodies are far too light to matter; "
         "the binary star, planet-with-moon and small black hole presets are heavy enough to show it.");
    ImGui::Checkbox("Paths in view", &s.showPaths);
    help("Draws each body's trail, and its predicted future (dashed), into the scene. They are ray "
         "traced like everything else, so a path behind the hole appears bent around it.");
    if (s.showPaths) {
        ImGui::SameLine();
        ImGui::Checkbox("Predicted", &s.showPredictions);
        if (s.showPredictions) {
            ImGui::SliderFloat("Predict ahead", &s.predictionTime, 50.0f, 3000.0f, "%.0f x", ImGuiSliderFlags_Logarithmic);
            help("How far ahead, in units of the scene's dynamical time (M for a black hole).");
        }
        ImGui::SliderFloat("Line width", &s.pathWidth, 0.5f, 4.0f, "%.1f px");
    }

    ImGui::SeparatorText("Objects");
    presetButtons(ctx, true);
    ImGui::SeparatorText("Orbits");
    presetButtons(ctx, false);

    ctx.state.plannerOpen = false;
    // Open by default around a black hole; around a star or planet its planned orbit
    // would be one more ring in the view, so it starts closed there.
    static int plannerHost = -1;
    const int host = static_cast<int>(central::hostOf(s.central));
    if (host != plannerHost) {
        ImGui::SetNextItemOpen(host == static_cast<int>(central::Host::BlackHole));
        plannerHost = host;
    }
    if (ImGui::TreeNodeEx("Launch planner", ImGuiTreeNodeFlags_DefaultOpen)) {
        plannerSection(ctx);
        ImGui::TreePop();
    }

    if (ctx.particles.active()) {
        ImGui::SeparatorText("Debris");
        ImGui::Text("%d particles", ctx.particles.alive());
        ImGui::SameLine();
        if (ImGui::SmallButton("Clear")) {
            ctx.particles.clear();
            sys.tidalEvents.clear();  // the flare and outflow go with the debris
        }
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
            help("How long gas at 10,000 K takes to radiate its heat away; hotter gas cools faster (~T^-3). In units "
                 "of the scene's dynamical time around a big object.");
            ImGui::SliderFloat("Coolest", &s.gasFloorTemperature, 1000.0f, 15000.0f, "%.0f K", ImGuiSliderFlags_Logarithmic);
            help("The light around the hole keeps the gas from cooling below this.");
            ImGui::SliderFloat("Glows as", &s.debrisGlowFloor, 0.0f, 20000.0f, s.debrisGlowFloor > 0.0f ? "%.0f K" : "its own heat");
            help("Display only: the gas glows at least as if it were this hot. A torn-apart star's stream is lit by its "
                 "leftover heat and then by the flare; its pressure still follows its own temperature.");
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
                if (s.central.hasHorizon()) ImGui::Text("%.4g", b.state.r);
                else ImGui::TextUnformatted(central::formatLength(b.state.r, s.central.massSolar).c_str());
                ImGui::TableNextColumn();
                ImGui::Text("%.6f", b.constants.E);
                ImGui::SetItemTooltip("Energy per unit rest mass, conserved. Below 1 means bound.");
                ImGui::TableNextColumn();
                ImGui::Text("%.4g", b.constants.L);
                ImGui::SetItemTooltip("Angular momentum about the spin axis, conserved.");
                ImGui::TableNextColumn();
                ImGui::Text("%.4g", b.constants.Q);
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

// --- Tidal disruption flare -----------------------------------------------------------

// The light curve of a tidal disruption on log-log axes: the fallback rate (orange), the
// luminosity that escapes (blue; photon trapping flattens it near the peak), the
// t^(-5/3) law through the peak (grey, dashed), and now (white).
void lightCurve(const tde::Aftermath& after, double now) {
    const tde::Event& e = after.reference();
    const double t0 = e.time;
    const double tLo = std::max(0.5 * e.fallback.firstReturn, 1.0);
    const double tHi = std::max({8.0 * (after.peakTime() - t0), 1.5 * (now - t0), 4.0 * tLo});
    const double yLo = -3.0, yHi = 0.25;  // log10 of rate / peak

    const float em = ImGui::GetFontSize();
    ImVec2 origin = ImGui::GetCursorScreenPos();
    ImVec2 size(ImGui::GetContentRegionAvail().x, 8.0f * em);
    ImGui::InvisibleButton("lightcurve", size);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 corner(origin.x + size.x, origin.y + size.y);
    dl->AddRectFilled(origin, corner, IM_COL32(10, 11, 16, 255));
    dl->PushClipRect(origin, corner, true);
    auto px = [&](double t) { return origin.x + static_cast<float>(std::log(t / tLo) / std::log(tHi / tLo)) * size.x; };
    auto py = [&](double v) {
        const double y = std::log10(std::max(v, 1e-12));
        return corner.y - static_cast<float>((y - yLo) / (yHi - yLo)) * size.y;
    };
    // Decades of time.
    for (double d = std::pow(10.0, std::ceil(std::log10(tLo))); d < tHi; d *= 10.0) {
        dl->AddLine(ImVec2(px(d), origin.y), ImVec2(px(d), corner.y), IM_COL32(255, 255, 255, 22));
    }
    const int n = 220;
    std::vector<ImVec2> rate, light;
    for (int i = 0; i < n; ++i) {
        const double t = tLo * std::pow(tHi / tLo, static_cast<double>(i) / (n - 1));
        const double level = after.level(t0 + t);
        if (level > 1e-4) rate.emplace_back(px(t), py(level));
        const double lum = after.luminosity(t0 + t);
        if (lum > 1e-4) light.emplace_back(px(t), py(lum));
    }
    // The t^(-5/3) law through the peak.
    const double tp = after.peakTime() - t0;
    for (int i = 0; i < 24; ++i) {
        const double ta = tp * std::pow(1.25, 2 * i), tb = tp * std::pow(1.25, 2 * i + 1);
        dl->AddLine(ImVec2(px(ta), py(std::pow(ta / tp, -5.0 / 3.0))), ImVec2(px(tb), py(std::pow(tb / tp, -5.0 / 3.0))),
                    IM_COL32(150, 150, 160, 160), 1.0f);
    }
    if (light.size() > 1) dl->AddPolyline(light.data(), static_cast<int>(light.size()), IM_COL32(110, 170, 255, 255), 0, 1.5f);
    if (rate.size() > 1) dl->AddPolyline(rate.data(), static_cast<int>(rate.size()), IM_COL32(255, 160, 60, 255), 0, 2.0f);
    if (now > t0 + tLo) dl->AddLine(ImVec2(px(now - t0), origin.y), ImVec2(px(now - t0), corner.y), IM_COL32(255, 255, 255, 200));
    dl->PopClipRect();
    dl->AddText(ImVec2(origin.x + 4, origin.y + 2), IM_COL32(255, 160, 60, 255), "fallback");
    dl->AddText(ImVec2(origin.x + 4, origin.y + 2 + em), IM_COL32(110, 170, 255, 255), "light");
    dl->AddText(ImVec2(origin.x + 4, origin.y + 2 + 2 * em), IM_COL32(150, 150, 160, 255), "t^-5/3");
    char label[96];
    std::snprintf(label, sizeof(label), "time since disruption (log), to %.3g M", tHi);
    dl->AddText(ImVec2(corner.x - ImGui::CalcTextSize(label).x - 4, corner.y - em - 2), IM_COL32(150, 150, 160, 255), label);
}

void flareSection(UiContext& ctx) {
    Settings& s = ctx.settings;
    BodySystem& sys = ctx.bodies;
    ImGui::Checkbox("Flare and outflow", &s.tdeFlare);
    help("When a star is torn apart, the bound part of its debris swings out and comes back one orbit later, each "
         "parcel after its own exact (Kerr) radial period. Returning, the stream crosses itself (relativistic "
         "precession swings it round), shocks and settles into a small new disk: the hole flares. Much of the "
         "returning gas is blown out instead, by the shocks and the trapped radiation, in a wind that is wide while "
         "the fallback is strong and narrows into a funnel along the axis as it fades. This is what ESO's "
         "animation of AT2019qiz shows.");
    if (!sys.tidalEvents.empty()) {
        ImGui::SameLine();
        if (ImGui::SmallButton("Forget")) sys.tidalEvents.clear();
        ImGui::SetItemTooltip("Forget the disruptions so far (their debris stays).");
    }
    if (sys.tidalEvents.empty()) {
        ImGui::TextWrapped("No star has been torn apart yet. Try Famous black holes -> Tidal disruption, or "
                           "Objects -> Tidal disruption flare.");
    } else if (s.tdeFlare) {
        const tde::Aftermath after(sys.tidalEvents, outflowParams(s), s.central.massSolar, s.spin);
        if (after.active()) {
            const tde::Event& e = after.reference();
            const double now = sys.time;
            const double level = after.level(now);
            if (now < after.start()) {
                ImGui::Text("Debris coming back in %s", timeText(s, after.start() - now).c_str());
            } else {
                ImGui::Text("Fallback at %.0f%% of its peak, %s", 100.0 * level, now < after.peakTime() ? "rising" : "fading");
            }
            lightCurve(after, now);
            if (ImGui::BeginTable("flare", 2, ImGuiTableFlags_SizingStretchProp)) {
                auto row = [](const char* name, const std::string& value, const char* tip) {
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextDisabled("%s", name);
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(value.c_str());
                    if (tip) ImGui::SetItemTooltip("%s", tip);
                };
                row("Closest approach", lengthText(s, e.pericenter) + printf1(",  beta = %.2f", e.tidalRadius / e.pericenter),
                    "beta = r_t / r_p: how deep inside its tidal radius the star went.");
                row("Debris", printf1("%.0f%% bound, ", 100.0 * e.fallback.boundFraction) +
                                  printf1("%.0f%% escaping", 100.0 * e.unboundFraction),
                    "From each parcel's orbit when the star came apart: the unbound part leaves for good in a long "
                    "tail; the bound part returns.");
                row("First back", timeText(s, e.fallback.firstReturn) + "  after the disruption",
                    "The most bound debris returns after its own orbital period.");
                row("Peak", timeText(s, after.peakTime() - e.time) + printf1(",  then ~ t^%.2f", e.fallback.lateSlope),
                    "A flat spread of debris energies gives t^(-5/3) (Rees 1988; Phinney 1989).");
                row("Returned so far", printf1("%.1f%% of the star", 100.0 * after.returned(now)), nullptr);
                if (now > after.start()) {
                    row("Outflow front", lengthText(s, after.front(now)),
                        "The leading edge of the wind, coasting out from the circularization radius.");
                    row("Opening", printf1("+-%.0f degrees", after.opening(now) / kDeg),
                        "Half-opening angle of the gas launched now: wide at first, when the returning stream's "
                        "self-crossing shock throws it out every which way; a narrow funnel once a disk has formed.");
                }
                row("Flare (Sun-like star)",
                    printf1("%.3g x Eddington", after.eddingtonRatio(now)) +
                        printf1("  (accretion %.3g x)", after.accretionRatio(now)),
                    "What escapes when this hole eats a Sun-like star, at this point of the fallback: above Eddington "
                    "most of the light is trapped in the inflow and what escapes grows only as 1 + ln(L / L_Edd).");
                ImGui::EndTable();
            }
            const tde::Scales sun = tde::physicalScales(s.central.massSolar, 1.0, central::kSunRadiusKm);
            const double edd = central::eddingtonRate(s.central.massSolar) * 3.15576e7 / 1.98847e33;  // Suns a year
            ImGui::TextWrapped("A Sun-like star torn apart by a hole of %.3g Suns starts coming back after %.3g days and "
                               "falls back at up to %.2g Suns a year, %.0f times the Eddington rate. (The simulated star "
                               "is bigger, so that it can be seen and its debris returns in minutes, not weeks.)",
                               s.central.massSolar, sun.tMin / 86400.0, sun.peakRate, sun.peakRate / std::max(edd, 1e-30));
        }
    }

    ImGui::SeparatorText("Outflow");
    ImGui::SliderFloat("Blown out", &s.outflowFraction, 0.0f, 0.95f, "%.2f");
    help("Share of the returning debris that the wind carries away; the rest is accreted and powers the flare. "
         "Simulations find a large part, perhaps most, is blown out.");
    ImGui::SliderFloat("Speed", &s.outflowSpeed, 0.02f, 0.3f, "%.3f c", ImGuiSliderFlags_Logarithmic);
    help("Terminal speed along the axis (slower toward the rim). AT2019qiz's outflow ran at about 10,000 km/s "
         "(0.03 c); radio and X-ray outflows of other events reach 0.1-0.2 c.");
    ImGui::SliderFloat("Opening", &s.outflowOpeningDeg, 5.0f, 60.0f, "+-%.0f deg");
    help("Half-opening angle of the funnel once the fallback has faded.");
    ImGui::SliderFloat("Widening", &s.outflowWidening, 1.0f, 4.0f, "%.1f x");
    help("How much wider the outflow is at first, before a disk has formed to collimate it: its first, widest "
         "gas makes the cap; the funnel that follows makes the stem.");
    ImGui::SliderFloat("Thickness##outflow", &s.outflowThickness, 0.05f, 5.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
    help("Opacity of the wind. A real one is so opaque (electron scattering) that its photosphere lies far out, "
         "hiding the inside; here it is kept translucent enough to see.");
    ImGui::SliderFloat("Lit", &s.outflowBrightness, 0.05f, 10.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
    help("How brightly the flare lights the outflow.");
    ImGui::Checkbox("New disk", &s.flareDisk);
    help("With the accretion disk off (a quiet hole), the returning gas forms its own small disk, out to about twice "
         "the pericenter (where gas keeping the stream's angular momentum circles), shining with the flare.");
    ImGui::SameLine();
    ImGui::Checkbox("Film lighting", &s.flareCinematic);
    help("The flare's light on the outflow fades as R^-0.7 instead of the true R^-2, so that its far reaches show, "
         "as in artists' renditions. Off: true to the inverse-square law.");
    ImGui::Checkbox("Follow the event", &s.flareCamera);
    help("The camera's distance follows the event: close for the star's approach and its disruption, then back out "
         "with the stream and the outflow. Zooming yourself turns it off.");
    ImGui::SameLine();
    ImGui::Checkbox("Skip the wait", &s.flareFastForward);
    help("Once the stream has left, nothing much changes for hundreds of M until the bound debris comes back: time "
         "runs faster until then.");
    if (ctx.pace > 1.01) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.35f, 1.0f), "fast-forwarding x%.1f", ctx.pace);
    }
    ImGui::Checkbox("Strands in the debris", &s.debrisFilaments);
    help("The tide stretches whatever structure the star had into fine strands along the stream. Each parcel remembers "
         "where it sat in the star, across the stretch; the strands are drawn from that, so they move with the gas. "
         "They need finer steps through the gas, so they cost speed.");
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
    // The map is drawn in scene units: around a big star everything is divided by S.
    const double S = ctx.scale;
    const float inv = static_cast<float>(1.0 / S);
    const double surface = surfaceOf(s);

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
    // Project a Cartesian point in M (z = spin axis) onto the map.
    auto project = [&](float x, float y, float z) {
        float u = x * inv, v = (side ? z : y) * inv;
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
        const float mx = (mouse.x - center.x) / scale, my = -(mouse.y - center.y) / scale;  // scene units
        if (ImGui::IsItemActivated()) {
            auto [r, theta, phi] = kerr::fromCartesian(a, mx * S, my * S, 0.0);
            if (r > central::boundaryRadius(a, surface) + 0.3 * S) {
                ui.mapAiming = true;
                ui.launchR = static_cast<float>(r / S);
                ui.launchPolarDeg = 90.0f;
                ui.launchAzimuthDeg = static_cast<float>(std::fmod(phi / kDeg + 360.0, 360.0));
                ui.launchTiltDeg = 0.0f;
            }
        }
        if (ui.mapAiming && ImGui::IsItemActive()) {
            const double phi = ui.launchAzimuthDeg * kDeg;
            const double sx = ui.launchR * std::cos(phi), sy = ui.launchR * std::sin(phi);
            // 0.03 c of speed per scene unit of drag (scaled like orbital speeds).
            const double perUnit = 0.03 / std::sqrt(S);
            double vx = (mx - sx) * perUnit, vy = (my - sy) * perUnit;
            double speed = std::hypot(vx, vy);
            if (speed > 0.95) {
                vx *= 0.95 / speed;
                vy *= 0.95 / speed;
            }
            const double vr = vx * std::cos(phi) + vy * std::sin(phi);
            const double vphi = -vx * std::sin(phi) + vy * std::cos(phi);
            const double r = ui.launchR * S;
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
        double inner = s.diskInnerAtIsco ? kerr::iscoRadius(a, s.diskPrograde) : s.diskInner * S;
        if (surface > 0.0) inner = std::max(inner, surface);
        const double outer = std::max<double>(s.diskOuter * S, inner * 1.2);
        ImU32 diskCol = IM_COL32(255, 150, 60, 50);
        if (side) {
            float x0 = static_cast<float>(std::sqrt(inner * inner + a * a));
            float x1 = static_cast<float>(std::sqrt(outer * outer + a * a));
            dl->AddLine(project(x0, 0, 0), project(x1, 0, 0), IM_COL32(255, 150, 60, 140), 3.0f);
            dl->AddLine(project(-x0, 0, 0), project(-x1, 0, 0), IM_COL32(255, 150, 60, 140), 3.0f);
        } else {
            float rIn = static_cast<float>(std::sqrt(inner * inner + a * a)) * inv * scale;
            float rOut = static_cast<float>(std::sqrt(outer * outer + a * a)) * inv * scale;
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
    const double rh = central::boundaryRadius(a, surface);
    if (surface <= 0.0) {
        shell([&](double th) { return kerr::ergosphereRadius(a, th); }, IM_COL32(120, 90, 200, 200), 1.5f, false);
        shell([&](double) { return rh; }, IM_COL32(0, 0, 0, 255), 0.0f, true);
        shell([&](double) { return rh; }, IM_COL32(200, 200, 220, 255), 1.5f, false);
    } else {
        // The star or planet itself, in its own color.
        const glm::vec3 hue = s.central.kind == central::Kind::RockyPlanet ? glm::vec3(0.3f, 0.55f, 0.9f)
                            : s.central.kind == central::Kind::GasGiant    ? glm::vec3(0.85f, 0.7f, 0.5f)
                                                                            : blackbody::hue(central::displayTemperature(
                                                                                  s.central.surfaceTemperature));
        shell([&](double) { return rh; }, toColor(hue, 0.55f), 0.0f, true);
        shell([&](double) { return rh; }, toColor(hue, 1.0f), 1.5f, false);
    }
    if (!side) {
        auto circle = [&](double r, ImU32 col) {
            if (surface > 0.0 && r < surface) return;
            dl->AddCircle(center, static_cast<float>(std::sqrt(r * r + a * a)) * inv * scale, col, 96, 1.0f);
        };
        circle(kerr::iscoRadius(a, true), IM_COL32(80, 220, 120, 160));
        circle(kerr::iscoRadius(a, false), IM_COL32(220, 140, 80, 110));
        circle(kerr::photonOrbitRadius(a, true), IM_COL32(255, 255, 140, 120));
        // Rotation direction.
        float rr = static_cast<float>(rh) * inv * scale + 10.0f;
        dl->PathArcTo(center, rr, -0.5f, 0.9f, 16);
        dl->PathStroke(IM_COL32(200, 200, 220, 120), 0, 1.5f);
    } else {
        float tip = static_cast<float>(rh) * inv * scale + 16.0f;
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
        dl->AddCircleFilled(q, std::max(3.0f, b.radius * inv * scale), toColor(hue));
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
    if (ui.inspectorOpen && ui.inspectValid) {
        std::vector<ImVec2> pts;
        for (const auto& p : ui.inspect.path) pts.push_back(projectPoint(p));
        if (pts.size() > 1) dl->AddPolyline(pts.data(), static_cast<int>(pts.size()), kRayColor, 0, 1.5f);
    }

    if (nearest) {
        char label[96];
        std::snprintf(label, sizeof(label), "%s  r = %.4g M", nearest->name.c_str(), nearest->state.r);
        dl->AddText(ImVec2(nearestPos.x + 8, nearestPos.y - 8), IM_COL32(230, 230, 240, 255), label);
    }

    // Camera
    const Camera& cam = ctx.camera;
    auto cp = kerr::toCartesian(a, cam.r, cam.theta, cam.phi);
    ImVec2 c = project(static_cast<float>(cp[0]), static_cast<float>(cp[1]), static_cast<float>(cp[2]));
    dl->AddCircle(c, 5.0f, IM_COL32(120, 200, 255, 255), 12, 2.0f);
    dl->AddText(ImVec2(c.x + 7, c.y - 7), IM_COL32(120, 200, 255, 255), "camera");

    dl->PopClipRect();
    char caption[128];
    std::snprintf(caption, sizeof(caption), "%s%s",
                  side ? "side view: spin axis up" : ui.plannerOpen ? "top view: drag from a point to throw a body"
                                                                    : "top view: spin axis toward you",
                  S > 1.0 ? printf1("   (1 square unit = %.3g M)", S).c_str() : "");
    dl->AddText(ImVec2(origin.x + 6, maxCorner.y - ImGui::GetFontSize() - 4), IM_COL32(150, 150, 160, 255), caption);
    ImGui::End();
}

// --- Ray inspector -----------------------------------------------------------------

void inspectorWindow(UiContext& ctx) {
    UiState& ui = ctx.state;
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
            case trajectory::RayEnd::Body:
            case trajectory::RayEnd::Surface: {
                const char* what = r.end == trajectory::RayEnd::Disk ? "disk gas" : ui.inspectBody.c_str();
                ImGui::TextWrapped("Emitted by %s at r = %.4g M. Frequency shift g = %.4g: it arrives %s (%s).",
                                   what, r.rEnd, r.g, r.g >= 1.0 ? "bluer and brighter" : "redder and dimmer",
                                   r.g >= 1.0 ? "moving toward you outweighs gravity"
                                              : "gravity and/or receding motion");
                if (r.end == trajectory::RayEnd::Surface) {
                    ImGui::TextWrapped("It left %.1f degrees from the vertical, measured on the surface. Light bending "
                                       "lets light that leaves sideways (even backward-looking parts of the star) reach "
                                       "you.",
                                       std::acos(std::clamp(r.emissionCos, -1.0, 1.0)) * 180.0 / kerr::kPi);
                }
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
            case trajectory::RayEnd::Lens:
                ImGui::TextWrapped("Traced backwards it fell into %s: this pixel is part of its shadow.",
                                   ui.inspectBody.c_str());
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
        ImGui::Text("Went around the center %.2f times", r.turns);
        if (r.lensings > 0) ImGui::Text("Bent by %d small black hole%s", r.lensings, r.lensings > 1 ? "s" : "");
        ImGui::Text("Closest approach  r = %.4g M", r.rMin);
        ImGui::Text("Travel time       %s", timeText(ctx.settings, r.delay).c_str());
        ImGui::Text("Impact parameter  b = %.4g M", r.b);
        ImGui::Text("Carter constant   q = %.4g M^2", r.q);
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
            "Left drag      orbit the camera around the center\n"
            "Right drag     look around\n"
            "Scroll         move closer / farther\n"
            "W / S          move closer / farther\n"
            "A / D, Q / E   orbit sideways / up and down\n"
            "Ctrl+click     inspect the light reaching a pixel\n"
            "  (or middle-click)\n"
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
            "Every pixel follows a light ray backwards through curved spacetime: around a rotating (Kerr) "
            "black hole, or outside a neutron star, white dwarf, star or planet, which has the same "
            "spacetime outside its surface. For a black hole, the dark region is the shadow: rays that end "
            "up in the horizon. The thin bright ring around it is light that circled the hole before "
            "escaping. The disk on the far side appears above and below it, bent over the top by "
            "gravitational lensing.");
    }
    ImGui::End();
}

}  // namespace

bool plannedLaunch(const UiState& ui, double spin, double scale, double surface, kerr::Constants& c, kerr::State& s) {
    double v[3];
    if (!plannedVelocity(ui, spin, scale, surface, v)) return false;
    kerr::launchFromZamo(spin, ui.launchR * scale, ui.launchPolarDeg * kDeg, ui.launchAzimuthDeg * kDeg, v, 1.0, c, s);
    return true;
}

// A scale bar at the bottom of the view: how long a length is at the central object's
// distance (bent light magnifies things near a black hole or neutron star a little).
void scaleBar(UiContext& ctx) {
    const Camera& cam = ctx.camera;
    const Settings& s = ctx.settings;
    if (cam.insideHorizon(s.spin)) return;
    const ImVec2 screen = ImGui::GetIO().DisplaySize;
    const double pixelAngle = 2.0 * std::tan(0.5 * cam.fovDeg * kDeg) / std::max(screen.y, 1.0f);
    const double kmPerPixel = central::lengthKm(cam.r * pixelAngle, s.central.massSolar);
    if (!(kmPerPixel > 0.0)) return;
    // A round length (1, 2 or 5 times a power of ten) about 120 pixels long.
    const double target = 120.0 * kmPerPixel;
    const double decade = std::pow(10.0, std::floor(std::log10(target)));
    double km = decade;
    for (double step : {2.0, 5.0, 10.0}) {
        if (decade * step <= target) km = decade * step;
    }
    const float length = static_cast<float>(km / kmPerPixel);
    ImDrawList* dl = ImGui::GetForegroundDrawList();
    const float em = ImGui::GetFontSize();
    const ImVec2 a(screen.x - length - 2.0f * em, screen.y - 1.6f * em);
    const ImVec2 b(screen.x - 2.0f * em, a.y);
    const ImU32 color = IM_COL32(230, 230, 235, 220);
    dl->AddLine(a, b, color, 2.0f);
    dl->AddLine(ImVec2(a.x, a.y - 0.4f * em), ImVec2(a.x, a.y + 0.4f * em), color, 2.0f);
    dl->AddLine(ImVec2(b.x, b.y - 0.4f * em), ImVec2(b.x, b.y + 0.4f * em), color, 2.0f);
    const std::string label = central::formatLength(central::kmToM(km, s.central.massSolar), s.central.massSolar) +
                              (s.central.hasHorizon() ? "" : "  at the center's distance");
    const float width = ImGui::CalcTextSize(label.c_str()).x;
    dl->AddText(ImVec2(b.x - width, a.y - 1.5f * em), color, label.c_str());
}

void drawUi(UiContext& ctx) {
    UiState& ui = ctx.state;
    helpWindow(ui);
    inspectorWindow(ctx);
    scaleBar(ctx);
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
        if (ImGui::CollapsingHeader("Central object", ImGuiTreeNodeFlags_DefaultOpen)) centralSection(ctx);
        if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) cameraSection(ctx);
        if (ImGui::CollapsingHeader("Accretion disk and jets", ImGuiTreeNodeFlags_DefaultOpen)) diskSection(ctx);
        ui.plannerOpen = false;
        if (ImGui::CollapsingHeader("Objects and orbits", ImGuiTreeNodeFlags_DefaultOpen)) objectsSection(ctx);
        if (ctx.settings.central.hasHorizon()) {
            // Open by itself once a star has been torn apart.
            static size_t eventsSeen = 0;
            if (ctx.bodies.tidalEvents.size() > eventsSeen) ImGui::SetNextItemOpen(true);
            eventsSeen = ctx.bodies.tidalEvents.size();
            if (ImGui::CollapsingHeader("Tidal disruption flare")) flareSection(ctx);
        }
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
