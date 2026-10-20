#include "eht_view.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace {

constexpr double kPi = 3.14159265358979323846;
// The field shown: 160 uas across, as the EHT's published images (M around a hole of
// unknown distance).
constexpr double kFieldUas = 160.0;
constexpr double kFieldM = 16.0;
// For a hole at no particular distance, M is taken to subtend this much (for the blur).
constexpr double kNominalScale = 4.0;
// On the sky (north up, east left) M87*'s approaching jet points to position angle 288
// degrees (east of north). Its spin axis points away from us, so the jet coming toward us
// runs along minus the projected axis, minus beta: turned by 108 degrees counterclockwise
// it lies at 288. Sgr A*'s orientation is not known well: shown with its axis up.
constexpr double kM87Rotation = 108.0;

template <class... Args>
std::string format(const char* fmt, Args... args) {
    char text[512];
    std::snprintf(text, sizeof(text), fmt, args...);
    return text;
}

}  // namespace

const char* EhtView::targetName(int target) {
    switch (target) {
        case M87: return "M87*";
        case SgrA: return "Sgr A*";
        default: return "this black hole";
    }
}

EhtView::~EhtView() {
    shutdown();
}

void EhtView::stopWorker() {
    m_cancel = true;
    if (m_worker.joinable()) m_worker.join();
    m_cancel = false;
}

void EhtView::shutdown() {
    stopWorker();
    if (m_textures[0]) glDeleteTextures(2, m_textures);
    m_textures[0] = m_textures[1] = 0;
    m_started = false;
    m_hasShown = false;
}

void EhtView::work(Options o) {
    const auto start = std::chrono::steady_clock::now();
    const int threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency() - 1));
    eht::View view;
    view.a = o.spin;
    view.inclination = o.inclinationDeg * kPi / 180.0;
    view.frequency = o.frequencyGHz * 1e9;
    double thetaG = kNominalScale, target = 1.0;
    if (o.target == ThisHole) {
        view.massSolar = o.massSolar;
        view.halfWidth = kFieldM;
    } else {
        const eht::Target& t = eht::targets()[static_cast<size_t>(o.target)];
        view.massSolar = t.massSolar;
        thetaG = eht::angularScale(t.massSolar, t.distancePc);
        target = t.flux;
        view.halfWidth = 0.5 * kFieldUas / thetaG;
    }
    eht::Flow flow;
    flow.T0 = o.temperature;
    flow.beta = o.beta;
    flow.height = o.height;
    // The density that gives the observed flux at 230 GHz (other frequencies keep it).
    {
        eht::View at230 = view;
        at230.frequency = 230e9;
        flow.n0 = eht::calibrate(at230, flow, target, thetaG, threads, &m_cancel);
    }
    auto finish = [&](int size, bool preview) {
        view.size = size;
        Result r;
        r.options = o;
        r.model = eht::render(view, flow, threads, &m_cancel);
        if (m_cancel) return false;
        r.blurred = eht::blur(r.model, o.beamUas / thetaG);
        r.curve = eht::criticalCurve(view.a, view.inclination, 360);
        r.thetaG = thetaG;
        r.n0 = flow.n0;
        r.flux = eht::flux(r.model, thetaG);
        const double shadow = std::sqrt(27.0);
        r.ringUas = eht::ringDiameter(r.blurred, eht::curveCenter(r.curve), shadow) * thetaG;
        for (int k = 0; k < 72; ++k) r.shadowUas += 2.0 * eht::criticalRadius(r.curve, 2.0 * kPi * k / 72.0) * thetaG / 72.0;
        double brightest = 0.0;
        for (double i : r.model.intensity) brightest = std::max(brightest, i);
        r.brightest = eht::brightnessTemperature(brightest, view.frequency);
        r.preview = preview;
        r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::lock_guard<std::mutex> lock(m_mutex);
        m_pending = std::move(r);
        m_hasPending = true;
        return true;
    };
    if (o.size > 64 && !finish(64, true)) {
        m_running = false;
        return;
    }
    finish(o.size, false);
    m_running = false;
}

void EhtView::update(const Options& options) {
    if (!m_started || !(options == m_options)) {
        stopWorker();
        m_options = options;
        m_started = true;
        m_running = true;
        m_worker = std::thread(&EhtView::work, this, options);
    }
    if (!m_running && m_worker.joinable()) m_worker.join();
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_hasPending) {
        upload(m_pending);
        m_shown = std::move(m_pending);
        m_hasPending = false;
        m_hasShown = true;
    }
}

void EhtView::finish() {
    if (m_worker.joinable()) m_worker.join();
    m_running = false;
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_hasPending) {
        upload(m_pending);
        m_shown = std::move(m_pending);
        m_hasPending = false;
        m_hasShown = true;
    }
}

void EhtView::upload(const Result& r) {
    if (!m_textures[0]) {
        glCreateTextures(GL_TEXTURE_2D, 2, m_textures);
    }
    const eht::Image* images[2] = {&r.model, &r.blurred};
    for (int k = 0; k < 2; ++k) {
        const eht::Image& image = *images[k];
        // Linear brightness, as the EHT shows it: the model up to its 99.5th percentile (its
        // thin photon ring would otherwise leave the rest dark), the blurred image to its peak.
        std::vector<double> sorted(image.intensity);
        std::sort(sorted.begin(), sorted.end());
        const double top = sorted.empty() ? 1.0
                                          : (k == 0 ? sorted[static_cast<size_t>(0.995 * (sorted.size() - 1))]
                                                    : sorted.back());
        std::vector<unsigned char> rgba(image.intensity.size() * 4);
        for (size_t i = 0; i < image.intensity.size(); ++i) {
            float c[3];
            eht::afmhot(top > 0.0 ? image.intensity[i] / top : 0.0, c);
            for (int j = 0; j < 3; ++j) rgba[4 * i + j] = static_cast<unsigned char>(std::lround(255.0f * c[j]));
            rgba[4 * i + 3] = 255;
        }
        GLuint tex = m_textures[k];
        // A new size needs new storage: textures made with glTextureStorage2D are immutable.
        GLint width = 0;
        glGetTextureLevelParameteriv(tex, 0, GL_TEXTURE_WIDTH, &width);
        if (width != image.size) {
            glDeleteTextures(1, &tex);
            glCreateTextures(GL_TEXTURE_2D, 1, &tex);
            m_textures[k] = tex;
            glTextureStorage2D(tex, 1, GL_RGBA8, image.size, image.size);
            glTextureParameteri(tex, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTextureParameteri(tex, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTextureParameteri(tex, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTextureParameteri(tex, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        }
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTextureSubImage2D(tex, 0, 0, 0, image.size, image.size, GL_RGBA, GL_UNSIGNED_BYTE, rgba.data());
        glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    }
}

void EhtView::draw(float left, float bottomMargin) {
    const ImVec2 screen = ImGui::GetIO().DisplaySize;
    const float em = ImGui::GetFontSize();
    ImGui::SetNextWindowPos(ImVec2(left, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(std::max(screen.x - left, 10.0f), std::max(screen.y - bottomMargin, 10.0f)),
                             ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(1.0f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, IM_COL32(6, 6, 10, 255));
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                                   ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus |
                                   ImGuiWindowFlags_NoFocusOnAppearing;
    if (!ImGui::Begin("##telescope", nullptr, flags)) {
        ImGui::End();
        ImGui::PopStyleColor();
        return;
    }
    const Options& o = m_hasShown ? m_shown.options : m_options;
    ImGui::Text("What the Event Horizon Telescope sees: %s", targetName(o.target));
    ImGui::SameLine();
    ImGui::TextDisabled("  %.0f GHz,  spin a = %.2f,  seen %.0f deg from its axis%s", o.frequencyGHz, o.spin,
                        o.inclinationDeg, m_running ? "    (making the image...)" : "");
    if (!m_hasShown) {
        ImGui::TextDisabled("Ray tracing the flow's light...");
        ImGui::End();
        ImGui::PopStyleColor();
        return;
    }
    const Result& r = m_shown;
    const bool sky = o.target != ThisHole;
    const double rotation = (o.target == M87 ? kM87Rotation : 0.0) * kPi / 180.0;

    const ImVec2 origin = ImGui::GetCursorScreenPos();
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float gap = 1.5f * em;
    const float side = std::max(std::min((avail.x - gap) * 0.5f, avail.y - 7.5f * em), 4.0f * em);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const double hw = r.model.halfWidth;  // M
    const char* captions[2] = {"The model, ray traced", nullptr};
    char beamCaption[96];
    std::snprintf(beamCaption, sizeof(beamCaption), "As the EHT sees it: %.0f %s resolution", sky ? o.beamUas : o.beamUas / r.thetaG,
                  sky ? "uas" : "M");
    captions[1] = beamCaption;
    for (int k = 0; k < 2; ++k) {
        const ImVec2 corner(origin.x + k * (side + gap), origin.y);
        const ImVec2 center(corner.x + 0.5f * side, corner.y + 0.5f * side);
        const float scale = static_cast<float>(side / (2.0 * hw));  // pixels per M
        // Sky coordinates (M, beta up), turned by the position angle, to the screen.
        auto toScreen = [&](double alpha, double beta) {
            const double x = alpha * std::cos(rotation) - beta * std::sin(rotation);
            const double y = alpha * std::sin(rotation) + beta * std::cos(rotation);
            return ImVec2(center.x + static_cast<float>(x) * scale, center.y - static_cast<float>(y) * scale);
        };
        dl->PushClipRect(corner, ImVec2(corner.x + side, corner.y + side), true);
        dl->AddRectFilled(corner, ImVec2(corner.x + side, corner.y + side), IM_COL32(0, 0, 0, 255));
        dl->AddImageQuad(static_cast<ImTextureID>(m_textures[k]), toScreen(-hw, hw), toScreen(hw, hw), toScreen(hw, -hw),
                         toScreen(-hw, -hw), ImVec2(0, 0), ImVec2(1, 0), ImVec2(1, 1), ImVec2(0, 1));
        if (k == 0 && o.curve) {
            // The shadow's edge, the critical curve.
            std::vector<ImVec2> points;
            for (const eht::Point& p : r.curve) points.push_back(toScreen(p.alpha, p.beta));
            dl->AddPolyline(points.data(), static_cast<int>(points.size()), IM_COL32(90, 220, 255, 170),
                            ImDrawFlags_Closed, 1.2f);
        }
        if (k == 1) {
            // The beam, bottom right, as in the EHT's figures.
            const float beam = static_cast<float>(0.5 * o.beamUas / r.thetaG * scale);
            dl->AddCircle(ImVec2(corner.x + side - beam - 0.8f * em, corner.y + side - beam - 0.8f * em), beam,
                          IM_COL32(255, 255, 255, 200), 48, 1.5f);
            if (sky) {
                // The ring the EHT measured, dashed.
                const eht::Target& t = eht::targets()[static_cast<size_t>(o.target)];
                const float radius = static_cast<float>(0.5 * t.ringUas / r.thetaG * scale);
                for (int d = 0; d < 48; d += 2) {
                    const double a0 = 2.0 * kPi * d / 48.0, a1 = 2.0 * kPi * (d + 1) / 48.0;
                    dl->AddLine(ImVec2(center.x + radius * static_cast<float>(std::cos(a0)), center.y + radius * static_cast<float>(std::sin(a0))),
                                ImVec2(center.x + radius * static_cast<float>(std::cos(a1)), center.y + radius * static_cast<float>(std::sin(a1))),
                                IM_COL32(120, 255, 140, 200), 1.3f);
                }
            }
        }
        // Scale bar, bottom left: 40 uas, or 10 M.
        const double barM = sky ? 40.0 / r.thetaG : 10.0;
        const float bar = static_cast<float>(barM * scale);
        const ImVec2 b0(corner.x + 0.8f * em, corner.y + side - 1.0f * em);
        dl->AddLine(b0, ImVec2(b0.x + bar, b0.y), IM_COL32(255, 255, 255, 230), 2.0f);
        dl->AddText(ImVec2(b0.x, b0.y - 1.3f * em), IM_COL32(255, 255, 255, 230), sky ? "40 uas" : "10 M");
        if (sky && o.target == M87) {
            // Compass: north up, east left.
            const ImVec2 c0(corner.x + side - 1.2f * em, corner.y + 2.4f * em);
            dl->AddLine(c0, ImVec2(c0.x, c0.y - 1.4f * em), IM_COL32(255, 255, 255, 200), 1.5f);
            dl->AddLine(c0, ImVec2(c0.x - 1.4f * em, c0.y), IM_COL32(255, 255, 255, 200), 1.5f);
            dl->AddText(ImVec2(c0.x - 0.3f * em, c0.y - 2.5f * em), IM_COL32(255, 255, 255, 200), "N");
            dl->AddText(ImVec2(c0.x - 2.3f * em, c0.y - 0.6f * em), IM_COL32(255, 255, 255, 200), "E");
        }
        dl->PopClipRect();
        dl->AddText(ImVec2(corner.x, corner.y + side + 0.3f * em), IM_COL32(220, 220, 225, 255), captions[k]);
    }
    ImGui::Dummy(ImVec2(2.0f * side + gap, side + 1.6f * em));

    // The numbers.
    if (sky) {
        const eht::Target& t = eht::targets()[static_cast<size_t>(o.target)];
        ImGui::Text("%s", format("Flux %.2f Jy at %.0f GHz (density calibrated to the observed %.1f Jy at 230 GHz: "
                                 "n0 = %.2g cm^-3)", r.flux, o.frequencyGHz, t.flux, r.n0).c_str());
        if (o.frequencyGHz == 230.0 && std::abs(r.flux / t.flux - 1.0) > 0.1) {
            ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.3f, 1.0f),
                               "No density gives the observed flux at this electron temperature: the flow turns opaque. "
                               "Raise the temperature.");
        }
        ImGui::Text("%s", format("Ring %.1f uas across (the EHT measured %.0f +- %.1f); the shadow %.1f uas across; "
                                 "M subtends %.2f uas", r.ringUas, t.ringUas, t.ringErrorUas, r.shadowUas, r.thetaG).c_str());
        ImGui::Text("%s", format("Brightest pixel: %.2g K (brightness temperature)", r.brightest).c_str());
        ImGui::TextDisabled("The real image: %s, %s", t.reference, t.url);
    } else {
        ImGui::Text("%s", format("Ring %.1f M across; the shadow %.2f M (2 sqrt 27 = %.2f without spin). Brightest "
                                 "pixel %.2g K.", r.ringUas / r.thetaG, r.shadowUas / r.thetaG, 2.0 * std::sqrt(27.0),
                                 r.brightest).c_str());
    }
    char text[400];
    std::snprintf(text, sizeof(text),
                  "Thermal synchrotron light of a hot, thick flow (T_e = %.2g K (r/M)^-0.84, plasma beta %.0f) rotating "
                  "at 0.9 of Keplerian, falling in inside the ISCO. Cyan: the critical curve (the shadow's edge)%s. "
                  "%d x %d pixels%s, %.1f s.",
                  o.temperature, o.beta, sky ? "; green: the ring the EHT measured" : "", r.model.size, r.model.size,
                  r.preview ? " (a first look)" : "", r.seconds);
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s", text);
    ImGui::PopTextWrapPos();
    ImGui::End();
    ImGui::PopStyleColor();
}
