#include "merger_scene.h"

#include "film.h"

#include "blackbody.h"
#include "renderer.h"

#include <imgui.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace merger;

namespace {

double smooth01(double x) {
    x = std::clamp(x, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

// A round number of km near x, for grid lines and scale bars.
double niceKm(double x) {
    const double p = std::pow(10.0, std::floor(std::log10(std::max(x, 1e-9))));
    const double f = x / p;
    return (f < 1.5 ? 1.0 : f < 3.5 ? 2.0 : f < 7.5 ? 5.0 : 10.0) * p;
}

// Physical time as people say it: ms, s, minutes, hours, days.
std::string formatTime(double t) {
    char buf[64];
    const double a = std::abs(t);
    const char* sign = t < 0.0 ? "-" : "+";
    if (a < 1.0) std::snprintf(buf, sizeof(buf), "%s%.2f ms", sign, a * 1e3);
    else if (a < 120.0) std::snprintf(buf, sizeof(buf), "%s%.2f s", sign, a);
    else if (a < 7200.0) std::snprintf(buf, sizeof(buf), "%s%.1f minutes", sign, a / 60.0);
    else if (a < 2.0 * kDay) std::snprintf(buf, sizeof(buf), "%s%.1f hours", sign, a / 3600.0);
    else std::snprintf(buf, sizeof(buf), "%s%.1f days", sign, a / kDay);
    return buf;
}

std::string formatKm(double km) {
    char buf[64];
    if (km < 1e4) std::snprintf(buf, sizeof(buf), "%.0f km", km);
    else if (km < 1e8) std::snprintf(buf, sizeof(buf), "%.0f thousand km", km / 1e3);
    else std::snprintf(buf, sizeof(buf), "%.1f AU", km / 1.495978707e8);
    return buf;
}

}  // namespace

// --- Setup ----------------------------------------------------------------------------------------

bool MergerScene::init(const std::filesystem::path& shaderDir, std::string& error) {
    bool ok = m_sim.init(shaderDir, error);
    ok &= m_view.build(shaderDir, "fullscreen.vert", "merger_view.frag");
    ok &= m_map.build(shaderDir, "fullscreen.vert", "merger_map.frag");
    ok &= m_film.build(shaderDir, "fullscreen.vert", "merger_film.frag");
    ok &= m_noiseBake.buildCompute(shaderDir, "merger_noise.comp");
    if (!ok) error = errors();
    glCreateTextures(GL_TEXTURE_3D, 1, &m_noiseTexture);
    glTextureStorage3D(m_noiseTexture, 8, GL_RGBA16F, kNoiseDims, kNoiseDims, kNoiseDims);
    glTextureParameteri(m_noiseTexture, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTextureParameteri(m_noiseTexture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    for (GLenum wrap : {GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T, GL_TEXTURE_WRAP_R})
        glTextureParameteri(m_noiseTexture, wrap, GL_REPEAT);
    bakeNoise();
    glCreateTextures(GL_TEXTURE_2D, 1, &m_waveTexture);
    glTextureStorage2D(m_waveTexture, 1, GL_RG32F, kWaveSamples, 1);
    glTextureParameteri(m_waveTexture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(m_waveTexture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(m_waveTexture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(m_waveTexture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glCreateTextures(GL_TEXTURE_2D, 1, &m_skyTexture);
    glTextureStorage2D(m_skyTexture, 1, GL_R32F, 32, 16);
    glTextureParameteri(m_skyTexture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(m_skyTexture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(m_skyTexture, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTextureParameteri(m_skyTexture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    for (const EjectaComponent& c : gw170817Ejecta()) m_kilonova.emplace_back(c);
    return ok;
}

void MergerScene::shutdown() {
    m_audio.close();
    m_sim.shutdown();
    if (m_waveTexture) glDeleteTextures(1, &m_waveTexture);
    if (m_skyTexture) glDeleteTextures(1, &m_skyTexture);
    if (m_noiseTexture) glDeleteTextures(1, &m_noiseTexture);
    m_waveTexture = m_skyTexture = m_noiseTexture = 0;
}

void MergerScene::hotReload() {
    m_sim.hotReload();
    m_view.reloadIfChanged();
    m_map.reloadIfChanged();
    m_film.reloadIfChanged();
    if (m_noiseBake.reloadIfChanged()) bakeNoise();
}

std::string MergerScene::errors() const {
    return m_sim.errors() + m_view.error() + m_map.error() + m_film.error() + m_noiseBake.error();
}

void MergerScene::bakeNoise() {
    if (!m_noiseBake.valid() || !m_noiseTexture) return;
    m_noiseBake.use();
    glBindImageTexture(0, m_noiseTexture, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA16F);
    glDispatchCompute(kNoiseDims / 4, kNoiseDims / 4, kNoiseDims / 4);
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glGenerateTextureMipmap(m_noiseTexture);
}

void MergerScene::start(Camera& camera) {
    m_sim.reset(params);
    m_active = true;
    // The loudest the wave gets, about: two point masses a little over two radii apart
    // (2 mu M / a in d^2 I / dt^2). The simulation's own peak replaces it if larger.
    {
        const double M = params.m1 + params.m2, mu = params.m1 * params.m2 / M;
        m_wavePeak = 2.0 * mu * M / (2.2 * params.radiusKm / kKm);
    }
    m_skyUploaded = false;
    m_trace.clear();
    m_boomAge = 1e9;
    m_burstSeen = false;
    m_audioPhase = m_audioFreq = m_audioAmp = 0.0;
    m_wallTime = 0.0;
    m_orbitTime = 0.0;
    m_contactPhase = 0.0;
    setFilmTime(0.0);
    camera.reset();
    camera.scale = 1e15;  // the kilonova is drawn from billions of km away
    camera.fovDeg = 45.0f;
    camera.theta = 0.72;
    camera.phi = -1.2;
    camera.r = std::max(3.6 * m_sim.separation(), 8.0 * m_sim.radius());
    m_viewHalf = 0.75 * m_sim.separation() + 1.5 * m_sim.radius();
    m_boxHalf = 0.5 * m_sim.separation() + 1.8 * m_sim.radius();
    m_lastR = camera.r;
    m_lastTheta = camera.theta;
    m_lastPhi = camera.phi;
    m_autoAngles = true;
    autoCamera = true;
    const float black = 0.0f;
    glClearTexImage(m_skyTexture, 0, GL_RED, GL_FLOAT, &black);
}

void MergerScene::stop() {
    m_active = false;
    m_audio.close();
}

double MergerScene::sinceMerger() const {
    return m_sim.sinceMerger() * kSecond;
}

double MergerScene::burstTime() const {
    return breakoutTime(jet, kEjectaEdge * kBluePolarStretch * gw170817Ejecta()[0].velocity);
}

double MergerScene::ejectaFront(double t) const {
    return kEjectaEdge * kBluePolarStretch * gw170817Ejecta()[0].velocity * std::max(t, 0.0) / kSecond;
}

double MergerScene::kilonovaRamp(double t) const {
    // The winds that make most of the kilonova's mass blow off the hot remnant and the torus
    // over the first few hundred milliseconds.
    return smooth01((t - 0.006) / 0.15);
}

double MergerScene::cloudDepth(double t) const {
    if (!cloudOn || t <= 0.0) return 0.0;
    // The true cloud is opaque (debrisOpticalDepth: ~10^15 at first) and would hide the
    // remnant at once. It is drawn translucent instead, thinning as the true one does, as
    // 1 / R^2 (its column density), from 3 optical depths when it reaches 115 km (no more
    // before, so the collision shows through it), and let go once it is thinner than a few
    // tenths (by ~5 ms, 450 km out): later the kilonova's components show the ejecta. It
    // comes out over the first ~0.1 ms.
    const double R = debrisRadiusKm(debris, t);
    const double tau = std::min(4.0 * (100.0 / R) * (100.0 / R), 3.0);
    const double launch = smooth01(t / 1.2e-4);
    const double fade = smooth01((tau - 0.06) / 0.25);
    return tau * launch * fade;
}

double MergerScene::mergeFlash(double t) const {
    // The cores collide at a tenth of the speed of light and shock-heat to ~10^11 K: a flash
    // (mostly neutrinos in truth) rising within 40 microseconds and fading over ~0.4 ms.
    if (t <= 0.0 || !m_sim.merged()) return 0.0;
    return smooth01(t / 4e-5) * std::exp(-t / 4e-4);
}

double MergerScene::scale(const Camera& camera) const {
    return std::max(camera.r / 30.0, 1e-3);
}

// --- Time ---------------------------------------------------------------------------------------------

void MergerScene::setFilmTime(double t) {
    m_filmTime = std::fmod(std::max(t, 0.0), film::kLength);
    m_filmAudioTime = m_filmTime;
}

void MergerScene::update(double wallDt, bool paused, Camera& camera, int stepBudget, bool live) {
    if (!m_active) return;
    if (view == Film) {
        // The film keeps its own clock and does not need the hydrodynamics.
        if (!paused) m_filmTime = std::fmod(m_filmTime + wallDt * speed, film::kLength);
        synthesizeFilm(wallDt, live, paused);
        return;
    }
    if (!paused) m_wallTime += wallDt;
    if (!paused && sinceMerger() < kEndTime) {
        const double t = sinceMerger();
        const double collapse = params.collapseMs * 1e-3;
        const double rate = merger::tempo(tempo, t, collapse, jetOn ? burstTime() : -1.0) * speed;
        const double target = m_sim.time() + rate * wallDt / kSecond;
        m_sim.advance(target, stepBudget);
    }
    // Where the stars are as they touch: the tidal tails leave from the ends of that line.
    if (!m_sim.merged()) {
        const glm::dvec3 d = m_sim.starCenter(0) - m_sim.starCenter(1);
        if (d.x != 0.0 || d.y != 0.0) m_contactPhase = std::atan2(d.y, d.x);
    }
    // Once the hydrodynamics stops, where the ejecta went (for the kilonova's shape).
    if (!m_skyUploaded && m_sim.phase() == MergerSim::Phase::Coasting) {
        std::vector<float> sky = m_sim.ejectaSky();
        float peak = 0.0f;
        for (float f : sky) peak = std::max(peak, f);
        if (peak > 0.0f)
            for (float& f : sky) f /= peak;
        glTextureSubImage2D(m_skyTexture, 0, 0, 0, 32, 16, GL_RED, GL_FLOAT, sky.data());
        m_skyUploaded = true;
    }
    directCamera(wallDt, camera);
    synthesize(wallDt, live, paused);

    // The strain trace (face-on h+, normalized), one point per frame.
    const double now = m_sim.time();
    const auto I = m_sim.iddotAt(now);
    const double hp = 0.5 * (I[0] - I[1]);
    m_wavePeak = std::max(m_wavePeak, std::hypot(hp, I[3]));
    if (!paused) {
        m_trace.push_back(static_cast<float>(hp / m_wavePeak));
        if (m_trace.size() > 600) m_trace.erase(m_trace.begin());
    }
}

void MergerScene::directCamera(double wallDt, Camera& camera) {
    // You take over by zooming (distance) or dragging (angles).
    if (std::abs(camera.r / m_lastR - 1.0) > 1e-6) autoCamera = false;
    if (std::abs(camera.theta - m_lastTheta) > 1e-6 || std::abs(camera.phi - m_lastPhi) > 1e-6) m_autoAngles = false;
    const double t = sinceMerger();
    const double R = m_sim.radius();
    const bool merged = m_sim.merged();
    const bool coasting = m_sim.phase() == MergerSim::Phase::Coasting;

    // Where the action is (units): the stars, the remnant, then the ejecta.
    double simReach;
    if (!merged) simReach = 0.5 * m_sim.separation() + 1.6 * R;
    else simReach = std::clamp(0.9 * m_sim.extent() + 1.2 * R, 2.5 * R, 18.0 * R);
    // The fastest ejecta's front (the blue component's outer edge, as merger_view.frag draws
    // it), once the kilonova shows: the camera backs off to keep all of it in view.
    const double front = ejectaFront(t);
    const double shown = kilonovaOn ? kilonovaRamp(t) : 0.0;
    // While the debris cloud is thick the camera backs off to see it billow out around the
    // remnant, then closes in again on the remnant as it thins.
    const double cloudReach = debrisRadiusKm(debris, t) / kKm * smooth01((cloudDepth(t) - 0.06) / 0.4);
    // At 3.8 fronts the ejecta fill about a third of a 45-degree view, leaving room for the jet.
    const double target = std::max({3.3 * simReach, 2.8 * cloudReach, 3.8 * front * shown});
    if (autoCamera) {
        // Quick to back off (the ejecta grow exponentially on this clock), slow to close in.
        const double tau = target > camera.r ? 0.3 : 0.9;
        const double follow = 1.0 - std::exp(-wallDt / tau);
        camera.r = std::exp(std::log(camera.r) + (std::log(target) - std::log(camera.r)) * follow);
        // Never inside the cloud once it shows.
        if (shown > 0.03) camera.r = std::max(camera.r, 3.0 * front);
    }
    if (m_autoAngles) {
        // From above at first (the waves' spiral), then lower to see the torus and the jet,
        // rising and dipping as the camera circles so the scene is seen in depth: a full turn
        // every 45 s, up and down every 50 s.
        m_orbitTime += wallDt;
        const double lower = smooth01((t + 0.0005) / 0.004) * 0.35 + smooth01(std::log10(std::max(t, 1e-4) / 0.02)) * 0.2;
        const double rock = 0.2 * std::sin(2.0 * kPi * m_orbitTime / 50.0);
        const double thetaTarget = std::clamp(0.72 + lower + rock, 0.35, 1.48);
        camera.theta += (thetaTarget - camera.theta) * (1.0 - std::exp(-wallDt / 2.0));
        camera.phi += 2.0 * kPi / 45.0 * wallDt;
    }
    m_lastR = camera.r;
    m_lastTheta = camera.theta;
    m_lastPhi = camera.phi;

    // The drawing grid and the simulation view follow the matter.
    if (!coasting) {
        const double box = merged ? std::clamp(1.25 * m_sim.extent() + 2.0 * R, 2.5 * R, 30.0 * R)
                                  : 0.5 * m_sim.separation() + 1.8 * R;
        m_boxHalf += (box - m_boxHalf) * (1.0 - std::exp(-wallDt / 0.5));
    }
    double viewTarget = merged ? std::clamp(1.1 * m_sim.extent() + 1.5 * R, 3.0 * R, 40.0 * R)
                               : 0.75 * m_sim.separation() + 1.5 * R;
    if (coasting) viewTarget = std::max(viewTarget, 0.6 * (m_sim.time() - m_sim.collapseTime()) * 0.25 + 3.0 * R);
    viewTarget = std::min(viewTarget, 400.0 * R);
    m_viewHalf += (viewTarget - m_viewHalf) * (1.0 - std::exp(-wallDt / 0.8));
}

// --- Sound -------------------------------------------------------------------------------------------

void MergerScene::synthesize(double wallDt, bool live, bool paused) {
    m_lastAudio.clear();
    if (live && sound && !m_audio.isOpen()) m_audio.open(kSampleRate);
    if ((!live || !sound) && m_audio.isOpen()) m_audio.close();
    double want = wallDt * kSampleRate + m_sampleCarry;
    if (live) want += (0.12 - m_audio.queued()) * kSampleRate * 0.5;  // keep ~0.12 s queued
    const int n = std::max(0, static_cast<int>(want));
    m_sampleCarry = live ? 0.0 : want - n;
    if (n == 0) return;

    // The wave now: its frequency from how fast its phase turns, its amplitude relative to
    // the loudest so far.
    double freq = 0.0, amp = 0.0;
    if (!paused) {
        const double now = m_sim.time();
        const double dt = 2e-6 / kSecond;  // 2 microseconds
        const auto I1 = m_sim.iddotAt(now), I0 = m_sim.iddotAt(now - dt);
        const double p1 = std::atan2(I1[3], 0.5 * (I1[0] - I1[1])), p0 = std::atan2(I0[3], 0.5 * (I0[0] - I0[1]));
        const double f = std::abs(std::remainder(p1 - p0, 2.0 * kPi)) / (2.0 * kPi * dt * kSecond);  // Hz
        const double a = std::hypot(0.5 * (I1[0] - I1[1]), I1[3]) / m_wavePeak;
        freq = std::clamp(f * pitch, 25.0, 9000.0);
        amp = std::pow(std::clamp(a, 0.0, 1.0), 0.6) * volume * 0.5;
        // The hole has rung down: the waves are over.
        if (sinceMerger() > params.collapseMs * 1e-3 + 0.05) amp = 0.0;
    }
    // The gamma-ray burst: a deep rumble as the jet breaks out (a film's sound, not physics:
    // space carries no sound).
    const double t = sinceMerger();
    if (jetOn && !m_burstSeen && t > burstTime()) {
        m_burstSeen = true;
        m_boomAge = 0.0;
    }
    std::vector<int16_t> out(static_cast<size_t>(n));
    const double f0 = m_audioFreq > 0.0 ? m_audioFreq : freq, a0 = m_audioAmp;
    for (int k = 0; k < n; ++k) {
        const double u = (k + 1.0) / n;
        const double f = f0 + (freq - f0) * u;
        const double a = a0 + (amp - a0) * u;
        m_audioPhase = std::fmod(m_audioPhase + 2.0 * kPi * f / kSampleRate, 2.0 * kPi);
        // A little of the second harmonic, so low chirps are heard on small speakers.
        double s = a * (0.8 * std::sin(m_audioPhase) + 0.2 * std::sin(2.0 * m_audioPhase));
        if (m_boomAge < 6.0) {
            m_noise = m_noise * 1664525u + 1013904223u;
            const double white = (m_noise >> 9) / 4194304.0 - 1.0;
            m_lowpass += (white - m_lowpass) * 0.02;
            const double age = m_boomAge + static_cast<double>(k) / kSampleRate;
            s += volume * 6.0 * m_lowpass * std::exp(-age / 1.2) * std::min(age / 0.05, 1.0);
        }
        out[static_cast<size_t>(k)] = static_cast<int16_t>(std::clamp(s, -1.0, 1.0) * 32000.0);
    }
    m_boomAge += static_cast<double>(n) / kSampleRate;
    m_audioFreq = freq;
    m_audioAmp = amp;
    m_lastAudio = out;
    if (live && sound) m_audio.submit(out);
}

// The film's sound (a film's, not physics: space carries no sound): the gravitational-wave
// chirp, pitched up 40 times into hearing and following the stars' orbit exactly (twice its
// phase), rising in pitch and loudness to the merger; a boom at the flash, a rumble as the
// burst breaks out, rushes of air at the two whip pans, and a low hum under the calm last
// shot. Its clock is the film's, so it stays with the picture through pauses and loops.
void MergerScene::synthesizeFilm(double wallDt, bool live, bool paused) {
    m_lastAudio.clear();
    if (live && sound && !m_audio.isOpen()) m_audio.open(kSampleRate);
    if ((!live || !sound) && m_audio.isOpen()) m_audio.close();
    double want = wallDt * kSampleRate + m_sampleCarry;
    if (live) want += (0.12 - m_audio.queued()) * kSampleRate * 0.5;
    const int n = std::max(0, static_cast<int>(want));
    m_sampleCarry = live ? 0.0 : want - n;
    if (n == 0) return;
    // Film time across this buffer, from where the sound got to up to the picture.
    double t0 = m_filmAudioTime, t1 = m_filmTime;
    if (t1 < t0) t1 += film::kLength;  // looped
    if (paused) t1 = t0;
    std::vector<int16_t> out(static_cast<size_t>(n));
    constexpr double kChirpPitch = 40.0;
    double lastPsi = film::openingBinary(std::fmod(t0, film::kLength)).phase;
    for (int k = 0; k < n; ++k) {
        const double t = std::fmod(t0 + (t1 - t0) * (k + 1.0) / n, film::kLength);
        double s = 0.0;
        if (!paused) {
            m_noise = m_noise * 1664525u + 1013904223u;
            const double white = (m_noise >> 9) / 4194304.0 - 1.0;
            // The chirp, up to the flash.
            if (t < 6.62) {
                const double psi = film::openingBinary(t).phase;
                double dpsi = psi - lastPsi;
                if (dpsi < 0.0 || dpsi > 1.0) dpsi = 0.0;  // looped or jumped
                lastPsi = psi;
                m_filmChirpPhase = std::fmod(m_filmChirpPhase + 2.0 * kChirpPitch * dpsi, 2.0 * kPi);
                const double rate = std::min(dpsi * kSampleRate / std::max((t1 - t0) / n * kSampleRate, 1e-9), 40.0);
                const double loud = std::pow(std::clamp(rate / 30.0, 0.0, 1.0), 0.67) * (1.0 - smooth01((t - 6.3) / 0.32));
                s += 0.32 * loud * (0.8 * std::sin(m_filmChirpPhase) + 0.2 * std::sin(2.0 * m_filmChirpPhase));
            }
            // The boom at the flash and the rumble of the burst: low-passed noise and a thump.
            const double boom = t - 6.62, burst = t - 8.07;
            m_boomLow += (white - m_boomLow) * 0.015;
            m_rumbleLow += (white - m_rumbleLow) * 0.03;
            if (boom >= 0.0 && boom < 5.0)
                s += (5.0 * m_boomLow + 0.5 * std::sin(2.0 * kPi * 42.0 * boom)) * std::exp(-boom / 1.1) *
                     std::min(boom / 0.02, 1.0);
            if (burst >= 0.0 && burst < 4.0)
                s += 3.5 * m_rumbleLow * std::exp(-burst / 0.9) * std::min(burst / 0.05, 1.0);
            // Rushes of air at the whip pans.
            m_whooshLow += (white - m_whooshLow) * 0.25;
            for (double at : {9.97, 10.9}) {
                const double a = t - at + 0.15;
                if (a >= 0.0 && a < 0.45) s += 0.6 * (white - m_whooshLow) * std::sin(kPi * a / 0.45);
            }
            // The hum under the last shot.
            if (t > 14.9) {
                const double fade = smooth01((t - 14.9) / 0.8);
                s += 0.07 * fade * (std::sin(2.0 * kPi * 55.0 * t) + 0.6 * std::sin(2.0 * kPi * 82.5 * t));
            }
            s *= volume;
        }
        out[static_cast<size_t>(k)] = static_cast<int16_t>(std::clamp(std::tanh(s), -1.0, 1.0) * 32000.0);
    }
    m_filmAudioTime = std::fmod(t1, film::kLength);
    m_lastAudio = out;
    if (live && sound) m_audio.submit(out);
}

// --- Drawing ----------------------------------------------------------------------------------------------

void MergerScene::buildWaveTexture(double window) {
    // H(t) = (Iddot_xx - Iddot_yy) / 2 + i Iddot_xy, the wave going out in the orbital plane,
    // normalized to the loudest so far; oldest first.
    std::vector<float> data(2 * kWaveSamples);
    const double now = m_sim.time();
    for (int k = 0; k < kWaveSamples; ++k) {
        const double t = now - window * (1.0 - static_cast<double>(k) / (kWaveSamples - 1));
        const auto I = m_sim.iddotAt(t);
        data[2 * k] = static_cast<float>(0.5 * (I[0] - I[1]) / m_wavePeak);
        data[2 * k + 1] = static_cast<float>(I[3] / m_wavePeak);
    }
    glTextureSubImage2D(m_waveTexture, 0, 0, 0, kWaveSamples, 1, GL_RG, GL_FLOAT, data.data());
    m_waveWindow = window;
}

void MergerScene::render(Renderer& renderer, const Settings& settings, const Camera& camera, int width, int height,
                         float renderScale) {
    if (!m_active) return;
    if (view == Film) {
        renderer.beginExternal(width, height, renderScale);
        drawFilm(renderer, renderer.traceWidth(), renderer.traceHeight());
        // The film's own look: its shader writes the final colours (white-hot parts above 1
        // simply clip, as the reference's do), with a soft glow around everything bright.
        Settings post = settings;
        post.bloomStrength = 0.12f;
        post.exposure = 1.0f;
        post.tonemapper = 2;
        renderer.finishExternal(post, 1.0f);
        renderer.bindOutput();
        drawOverlays(width, height, width);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return;
    }
    // The simulation's grids around the action.
    const float mapHalf = static_cast<float>(m_viewHalf * 1.9);
    m_sim.deposit(glm::vec3(0.0f), static_cast<float>(m_boxHalf), glm::vec2(0.0f), mapHalf);

    // Ultra HD: the merger is traced at the screen's full resolution at least (unless the
    // resolution is left to follow a frame-time target).
    if (!settings.dynamicResolution) renderScale = std::max(renderScale, 1.0f);
    const int splitX = view == SideBySide ? width / 2 : (view == Cinematic ? width : 0);
    renderer.beginExternal(width, height, renderScale);
    if (view != Simulation) {
        const int th = renderer.traceHeight();
        const int tw = renderer.traceWidth();
        const int cw = view == SideBySide ? tw / 2 : tw;
        drawCinematic(renderer, settings, camera, 0, 0, cw, th);
    }
    const double t = sinceMerger();
    float flash = view != Simulation ? static_cast<float>(0.5 * mergeFlash(t)) : 0.0f;
    if (jetOn && t > burstTime()) flash += static_cast<float>(0.6 * std::exp(-(t - burstTime()) / 0.5));
    // A stronger glow than the black hole scenes: films show neutron stars haloed in light.
    Settings post = settings;
    post.bloomStrength = std::max(settings.bloomStrength, 0.14f);
    renderer.finishExternal(post, 1.0f + flash);
    renderer.bindOutput();
    if (view != Cinematic) drawMap(splitX, 0, width - splitX, height);
    drawOverlays(width, height, splitX);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void MergerScene::drawCinematic(Renderer& renderer, const Settings& settings, const Camera& camera, int x, int y, int w,
                                int h) {
    const double S = scale(camera);
    const double t = sinceMerger();
    const double R = m_sim.radius();
    ShaderProgram& p = m_view;
    p.use();
    glViewport(x, y, w, h);
    p.set("uResolution", glm::vec2(static_cast<float>(w), static_cast<float>(h)));
    p.set("uViewport", glm::vec2(static_cast<float>(x), static_cast<float>(y)));
    // The camera, in Cartesian render units.
    const double st = std::sin(camera.theta), ct = std::cos(camera.theta);
    const double sp = std::sin(camera.phi), cp = std::cos(camera.phi);
    const glm::vec3 er(static_cast<float>(st * cp), static_cast<float>(st * sp), static_cast<float>(ct));
    const glm::vec3 eth(static_cast<float>(ct * cp), static_cast<float>(ct * sp), static_cast<float>(-st));
    const glm::vec3 eph(static_cast<float>(-sp), static_cast<float>(cp), 0.0f);
    const glm::mat3 frame(er, eth, eph);
    p.set("uCamPos", er * static_cast<float>(camera.r / S));
    p.set("uCamBasis", frame * camera.basis());
    const float tanHalf = std::tan(0.5f * glm::radians(camera.fovDeg));
    p.set("uTanHalfFov", tanHalf);
    p.set("uPixelAngle", 2.0f * tanHalf / static_cast<float>(h));
    p.set("uTime", static_cast<float>(m_wallTime));
    p.set("uFilm", filmColor ? 1 : 0);

    // Lenses: the two stars, or the black hole.
    const float inv = static_cast<float>(1.0 / S);
    if (m_sim.hasHole()) {
        const merger::BlackHole& bh = m_sim.hole();
        p.set("uLens[0]", glm::vec4(glm::vec3(bh.position) * inv, static_cast<float>(m_sim.holeMass()) * inv));
        p.set("uLensCount", 1);
        p.set("uHole", 0);
    } else {
        for (int k = 0; k < 2; ++k) {
            const std::string name = "uLens[" + std::to_string(k) + "]";
            p.set(name.c_str(), glm::vec4(glm::vec3(m_sim.starCenter(k)) * inv, static_cast<float>(m_sim.starMass(k)) * inv));
        }
        p.set("uLensCount", 2);
        p.set("uHole", -1);
    }

    // The simulation's matter.
    const bool volumeOn = m_sim.boxSize() * inv > 1e-3f;
    p.set("uVolumeOn", volumeOn ? 1 : 0);
    p.set("uBoxMin", m_sim.boxMin() * inv);
    p.set("uBoxSize", m_sim.boxSize() * inv);
    p.set("uOpacityScale", static_cast<float>(12.0 * S));
    p.set("uMatterGlow", 1.0f);
    p.set("uLabelScale", static_cast<float>(S / R));
    glBindTextureUnit(1, m_sim.volumeA());
    glBindTextureUnit(2, m_sim.volumeB());

    // Gravitational waves in the orbital plane, out to near the camera.
    const bool wavesVisible = waves && t < params.collapseMs * 1e-3 + 0.01;
    p.set("uWavesOn", wavesVisible ? 1 : 0);
    if (wavesVisible) {
        const double outer = 0.9 * camera.r;
        buildWaveTexture(outer);
        p.set("uWaveNow", 1.0f);
        p.set("uWaveStep", static_cast<float>(S / outer));
        p.set("uWaveInner", static_cast<float>(1.3 * R / S));
        p.set("uWaveOuter", static_cast<float>(outer / S));
        p.set("uWaveGlow", 1.4f);
        glBindTextureUnit(3, m_waveTexture);
    }

    // The debris cloud: its edge, optical depth (as drawn), brightness and temperature (as
    // shown); where the remnant ends, the line the stars touched along, how far its billows
    // have churned, how brightly the remnant lights it from inside.
    const double depth = cloudDepth(t);
    glm::vec4 cloud(0.0f), cloudMore(0.0f);
    if (depth > 1e-3) {
        const double T = debrisTemperature(debris, t);
        // Its true temperature (10^9 - 10^11 K) is shown compressed, as 2,500 K at 10^9 K
        // and rising as its square root, so that it cools from blue-white through yellow and
        // orange to red as it expands.
        const double shownT = std::clamp(2500.0 * std::sqrt(T / 1e9), 1500.0, 20000.0);
        const double bright = std::clamp(std::pow(shownT / 9000.0, 0.8), 0.1, 1.6);
        cloud = glm::vec4(static_cast<float>(debrisRadiusKm(debris, t) / kKm / S), static_cast<float>(depth),
                          static_cast<float>(bright), static_cast<float>(shownT));
        const double churn = 0.6 * std::log1p(t / 3e-4) + 0.03 * m_wallTime;
        cloudMore = glm::vec4(static_cast<float>(1.1 * R / S), static_cast<float>(m_contactPhase),
                              static_cast<float>(churn), static_cast<float>(1.0 + 1.5 * mergeFlash(t)));
    }
    p.set("uCloud", cloud);
    p.set("uCloudMore", cloudMore);
    p.set("uMergeFlash", static_cast<float>(mergeFlash(t)));

    // The kilonova: each component's size, temperature, brightness and optical depth now.
    // (Not drawn until it shows: marching an invisible kilonova costs as much as a visible one.)
    const double ramp = kilonovaOn ? kilonovaRamp(t) : 0.0;
    const bool kilonovaShows = ramp > 0.003 && t > 0.0;
    p.set("uKilonovaOn", kilonovaShows ? 1 : 0);
    const auto ej = gw170817Ejecta();
    for (int k = 0; k < 3; ++k) {
        const std::string name = "uKilonova[" + std::to_string(k) + "]";
        glm::vec4 kn(0.0f);
        if (kilonovaShows) {
            const Kilonova& model = m_kilonova[static_cast<size_t>(k)];
            const double outer = kEjectaEdge * ej[k].velocity * t / kSecond / S;
            // What the eye sees is surface brightness, sigma T^4 at the photosphere (not the
            // total luminosity: the young ejecta are small but hot), compressed to fit.
            const double Tphot = model.temperature(t);
            const double T = std::clamp(Tphot, 1500.0, 30000.0);
            const double rCm = ej[k].velocity * kLightCm * t;
            const double tau = 3.0 * ej[k].opacity * ej[k].mass * kSunGrams / (4.0 * kPi * rCm * rCm);
            const double bright = 0.85 * std::clamp(Tphot / 9000.0, 0.03, 1.6);
            kn = glm::vec4(static_cast<float>(outer), static_cast<float>(T), static_cast<float>(bright * ramp),
                           static_cast<float>(std::min(tau, 60.0) * ramp));
        }
        p.set(name.c_str(), kn);
    }
    glBindTextureUnit(4, m_skyTexture);
    glBindTextureUnit(5, m_noiseTexture);

    // The jet and its cocoon.
    glm::vec4 jetU(0.0f);
    glm::vec2 jetMore(0.0f);
    float flash = 0.0f;
    if (jetOn && t > jet.launch) {
        const double tb = burstTime();
        const double c = 1.0 / kSecond / S;  // render units per second at the speed of light
        const double zBreak = jet.headSpeed * (tb - jet.launch) * c;
        const double head = t < tb ? jet.headSpeed * (t - jet.launch) * c : zBreak + (t - tb) * c;
        double bright = 0.5;
        if (t > tb) bright = 6.0 * std::exp(-(t - tb) / 1.0) + 1.2 / (1.0 + (t - tb) / 3.0);
        const double cocoon = t > tb ? jet.cocoonSpeed * (t - tb) * c : 0.0;
        jetU = glm::vec4(static_cast<float>(std::min(head, 4.0 * camera.r / S)), static_cast<float>(cocoon),
                         static_cast<float>(bright), static_cast<float>(1.3 * jet.coreAngle));
        jetMore = glm::vec2(static_cast<float>(zBreak), static_cast<float>(t > tb ? 2.2 / (1.0 + (t - tb) / 4.0) : 0.0));
        if (t > tb) flash = static_cast<float>(1.5 * std::exp(-(t - tb) / 0.7));
    }
    p.set("uJet", jetU);
    p.set("uJetMore", jetMore);
    p.set("uJetAxis", m_sim.hasHole() ? glm::vec3(m_sim.hole().axis) : glm::vec3(0.0f, 0.0f, 1.0f));
    p.set("uFlash", flash);

    glBindTextureUnit(0, renderer.blackbodyTexture());
    p.set("uBlackbodyRange", glm::vec2(static_cast<float>(blackbody::kLogTMin), static_cast<float>(blackbody::kLogTMax)));
    p.set("uStarBrightness", settings.starBrightness);
    p.set("uGalaxyBrightness", settings.galaxyBrightness);
    renderer.drawFullscreen();
}

void MergerScene::drawFilm(Renderer& renderer, int w, int h) {
    const double t = m_filmTime;
    const film::Binary b = t >= 14.8 ? film::closingBinary(t) : film::openingBinary(t);
    ShaderProgram& p = m_film;
    p.use();
    glViewport(0, 0, w, h);
    p.set("uResolution", glm::vec2(static_cast<float>(w), static_cast<float>(h)));
    // The 9:16 frame fills the height (or the width, if the view is narrower than it).
    p.set("uFrameHeight", static_cast<float>(std::min(static_cast<double>(h), w / film::kAspect)));
    p.set("uTime", static_cast<float>(t));
    p.set("uStar0", glm::vec3(film::starPosition(b, 0)));
    p.set("uStar1", glm::vec3(film::starPosition(b, 1)));
    p.set("uGlare", static_cast<float>(b.radius));
    p.set("uMerge", static_cast<float>(b.merge));
    p.set("uOrbitCentre", glm::vec3(b.centre));
    p.set("uOrbitPhase", static_cast<float>(b.phase));
    p.set("uOrbitRadius", static_cast<float>(b.orbit));
    renderer.drawFullscreen();
}

void MergerScene::drawMap(int x, int y, int w, int h) {
    ShaderProgram& p = m_map;
    p.use();
    glViewport(x, y, w, h);
    p.set("uViewport", glm::vec2(static_cast<float>(x), static_cast<float>(y)));
    p.set("uResolution", glm::vec2(static_cast<float>(w), static_cast<float>(h)));
    p.set("uCenter", glm::vec2(0.0f));
    // Once the hydrodynamics stops: edge-on, the ejecta.
    const double t = sinceMerger();
    const bool edgeOn = m_sim.phase() == MergerSim::Phase::Coasting && kilonovaOn;
    p.set("uEdgeOn", edgeOn ? 1 : 0);
    if (edgeOn) {
        const auto ej = gw170817Ejecta();
        double peak = -40.0;
        for (int k = 0; k < 3; ++k) {
            const std::string name = "uEjecta[" + std::to_string(k) + "]";
            const double R = kEjectaEdge * ej[k].velocity * t / kSecond;
            const double rCm = R * kKm * 1e5;
            const double mean = std::log10(ej[k].mass * kSunGrams / (4.0 / 3.0 * kPi * rCm * rCm * rCm));
            peak = std::max(peak, mean + 0.3);
            p.set(name.c_str(), glm::vec4(static_cast<float>(R), static_cast<float>(mean), 0.0f, 0.0f));
        }
        glm::vec4 funnel(0.0f);
        if (jetOn && t > jet.launch) {
            const double tb = burstTime();
            const double c = 1.0 / kSecond;
            const double zBreak = jet.headSpeed * (tb - jet.launch) * c;
            const double head = t < tb ? jet.headSpeed * (t - jet.launch) * c : zBreak + (t - tb) * c;
            funnel = glm::vec4(static_cast<float>(head), static_cast<float>(jet.coreAngle * 1.3),
                               static_cast<float>(t > tb ? jet.cocoonSpeed * (t - tb) * c : 0.0), static_cast<float>(zBreak));
        }
        p.set("uFunnel", funnel);
        p.set("uHalfHeight", static_cast<float>(1.15 * ejectaFront(t)));
        p.set("uLogRange", glm::vec2(static_cast<float>(peak - 4.5), static_cast<float>(peak + 0.4)));
        p.set("uGridStep", static_cast<float>(niceKm(ejectaFront(t) * kKm / 2.0) / kKm));
        glDrawArrays(GL_TRIANGLES, 0, 3);
        return;
    }
    p.set("uHalfHeight", static_cast<float>(m_viewHalf));
    p.set("uMapMin", m_sim.mapMin());
    p.set("uMapSize", m_sim.mapSize());
    // From the stars' centres (density x diameter) down six and a half decades.
    const double peak = std::log10(std::max(polytropeCentralDensity(params.m1, m_sim.radius()) * kGramsPerCc, 1.0) *
                                   1.3 * m_sim.radius() * kKm * 1e5);
    p.set("uLogRange", glm::vec2(static_cast<float>(peak - 6.0), static_cast<float>(peak + 0.6)));
    p.set("uGridStep", static_cast<float>(niceKm(m_viewHalf * kKm / 2.0) / kKm));
    glBindTextureUnit(0, m_sim.map());
    glDrawArrays(GL_TRIANGLES, 0, 3);
}

// --- Overlays -----------------------------------------------------------------------------------------------

void MergerScene::drawOverlays(int width, int height, int splitX) {
    // ImGui registers its font texture during its first frame's Render: until then there is
    // nothing to draw text with.
    if (ImGui::GetFrameCount() <= 1) return;
    // Built twice: text at a size not drawn before bakes new glyphs, which can grow the font
    // atlas and leave the first list's texture coordinates stale.
    {
        ImDrawList warm(ImGui::GetDrawListSharedData());
        warm._ResetForNewFrame();
        warm.PushClipRect(ImVec2(0.0f, 0.0f), ImVec2(static_cast<float>(width), static_cast<float>(height)));
        warm.PushTexture(ImGui::GetIO().Fonts->TexRef);
        fillOverlays(warm, width, height, splitX);
    }
    ImDrawList list(ImGui::GetDrawListSharedData());
    list._ResetForNewFrame();
    list.PushClipRect(ImVec2(0.0f, 0.0f), ImVec2(static_cast<float>(width), static_cast<float>(height)));
    list.PushTexture(ImGui::GetIO().Fonts->TexRef);
    fillOverlays(list, width, height, splitX);
    list.PopTexture();
    list.PopClipRect();
    ImDrawData dd;
    dd.Valid = true;
    dd.CmdLists.push_back(&list);
    dd.CmdListsCount = 1;
    dd.TotalVtxCount = list.VtxBuffer.Size;
    dd.TotalIdxCount = list.IdxBuffer.Size;
    dd.DisplayPos = ImVec2(0.0f, 0.0f);
    dd.DisplaySize = ImVec2(static_cast<float>(width), static_cast<float>(height));
    dd.FramebufferScale = ImVec2(1.0f, 1.0f);
    dd.Textures = &ImGui::GetPlatformIO().Textures;
    ImGui_ImplOpenGL3_RenderDrawData(&dd);
}

void MergerScene::fillOverlays(ImDrawList& list, int width, int height, int splitX) {
    if (view == Film) {
        if (captions) fillFilmCaptions(list, width, height);
        return;
    }
    ImFont* font = ImGui::GetFont();
    const float u = static_cast<float>(height) / 1080.0f;
    auto text = [&](float size, ImVec2 pos, ImU32 col, const std::string& s, float align = 0.0f) {
        const ImVec2 dim = font->CalcTextSizeA(size, 1e9f, 0.0f, s.c_str());
        pos.x -= dim.x * align;
        list.AddText(font, size, ImVec2(pos.x + 2.0f * u, pos.y + 2.0f * u), IM_COL32(0, 0, 0, 200), s.c_str());
        list.AddText(font, size, pos, col, s.c_str());
    };
    const double t = sinceMerger();
    const double tb = burstTime();

    // --- The cinematic side: a caption for what is happening, and the clock ---------------------
    if (view != Simulation && captions) {
        const float cx = 0.5f * static_cast<float>(splitX);
        std::string title, line;
        char buf[256];
        if (!m_sim.merged() && t < -0.0015) {
            title = "Two neutron stars spiral together";
            std::snprintf(buf, sizeof(buf), "%.2f and %.2f Suns, each %.0f km across, %s apart, circling %.0f times a second",
                          params.m1, params.m2, 2.0 * params.radiusKm, formatKm(m_sim.separation() * kKm).c_str(),
                          std::sqrt(m_sim.totalMass() / std::pow(m_sim.separation(), 3)) / (2.0 * kPi) / kSecond);
            line = buf;
        } else if (!m_sim.hasHole()) {
            if (t < 0.0012) {
                title = "They merge";
                line = "A flash as the cores collide; tides tear their outer layers into spiral arms";
            } else if (cloudDepth(t) > 0.3) {
                title = "A cloud of debris bursts out";
                std::snprintf(buf, sizeof(buf),
                              "%.3f Suns of neutron-rich gas flung out at up to %.1f c, cooling as it expands: %.1f billion K",
                              debris.mass, debris.edgeSpeed, debrisTemperature(debris, t) / 1e9);
                line = buf;
            } else {
                title = "A hypermassive neutron star";
                std::snprintf(buf, sizeof(buf), "%.2f Suns, too heavy to last, held up by its spin and heat: %.0f billion K",
                              m_sim.totalMass(), m_sim.maxTemperature() / 1e9);
                line = buf;
            }
        } else if (t < jet.launch) {
            title = "It collapses into a black hole";
            std::snprintf(buf, sizeof(buf), "%.2f Suns, spinning at %.2f of the maximum, inside a torus of %.2f Suns",
                          m_sim.holeMass(), m_sim.hole().spin, m_sim.hole().diskMass);
            line = buf;
        } else if (jetOn && t < tb) {
            title = "A jet drills out through the debris";
            line = "Launched along the black hole's spin axis at a large fraction of the speed of light";
        } else if (jetOn && t < tb + 6.0) {
            title = "Gamma-ray burst";
            std::snprintf(buf, sizeof(buf),
                          "The jet breaks out %.1f s after the merger (GRB 170817A came 1.7 s after GW170817)", tb);
            line = buf;
        } else if (t < 3600.0) {
            title = "Ejecting matter";
            double ejected = 0.0;
            for (const EjectaComponent& c : gw170817Ejecta()) ejected += c.mass;
            std::snprintf(buf, sizeof(buf), "%.2f Suns of neutron-rich debris flies out at up to a third of the speed of light",
                          ejected);
            line = buf;
        } else {
            title = "Kilonova: heavy elements glow";
            const Kilonova& blue = m_kilonova[0];
            const Kilonova& red = m_kilonova[2];
            std::snprintf(buf, sizeof(buf), "Gold, platinum and uranium forged and decaying: %.1e erg/s; blue side %.0f K, red %.0f K",
                          blue.luminosity(t) + m_kilonova[1].luminosity(t) + red.luminosity(t), blue.temperature(t),
                          red.temperature(t));
            line = buf;
        }
        text(46.0f * u, ImVec2(cx, 40.0f * u), IM_COL32(255, 255, 255, 255), title, 0.5f);
        text(24.0f * u, ImVec2(cx, 96.0f * u), IM_COL32(205, 225, 255, 235), line, 0.5f);
        text(30.0f * u, ImVec2(30.0f * u, static_cast<float>(height) - 70.0f * u), IM_COL32(255, 255, 255, 230),
             "t = " + formatTime(t) + (m_sim.merged() ? "" : "  (to the merger)"));
    }

    // --- The simulation side: the wave, the clock, the scale ----------------------------------------------
    if (view != Cinematic) {
        const float x0 = static_cast<float>(splitX);
        const float w = static_cast<float>(width - splitX);
        const float H = static_cast<float>(height);
        const bool edgeOn = m_sim.phase() == MergerSim::Phase::Coasting && kilonovaOn;
        char buf[160];
        if (!edgeOn) {
            // The gravitational wave, running down the left edge (newest at the bottom).
            const float gx = x0 + 70.0f * u, gTop = 150.0f * u, gBottom = H - 150.0f * u;
            text(20.0f * u, ImVec2(gx, gTop - 50.0f * u), IM_COL32(230, 150, 110, 255), "Gravitational", 0.5f);
            text(20.0f * u, ImVec2(gx, gTop - 28.0f * u), IM_COL32(230, 150, 110, 255), "waves", 0.5f);
            list.AddLine(ImVec2(gx, gTop), ImVec2(gx, gBottom), IM_COL32(150, 90, 70, 180), 1.5f * u);
            if (m_trace.size() > 2) {
                std::vector<ImVec2> pts;
                const size_t n = m_trace.size();
                for (size_t k = 0; k < n; ++k) {
                    const float y = gTop + (gBottom - gTop) * static_cast<float>(k) / static_cast<float>(599);
                    pts.push_back(ImVec2(gx + 45.0f * u * m_trace[k], y));
                }
                list.AddPolyline(pts.data(), static_cast<int>(pts.size()), IM_COL32(255, 140, 90, 255), 0, 2.2f * u);
            }
            // The clock: milliseconds since the run began, as simulation movies show.
            text(18.0f * u, ImVec2(x0 + w - 30.0f * u, 30.0f * u), IM_COL32(200, 200, 210, 230), "Time in milliseconds", 1.0f);
            std::snprintf(buf, sizeof(buf), "%.2f", m_sim.time() * kSecond * 1e3);
            text(34.0f * u, ImVec2(x0 + w - 30.0f * u, 54.0f * u), IM_COL32(255, 255, 255, 255), buf, 1.0f);
            std::snprintf(buf, sizeof(buf), "%s merger", formatTime(t).c_str());
            text(18.0f * u, ImVec2(x0 + w - 30.0f * u, 96.0f * u), IM_COL32(200, 200, 210, 230), buf, 1.0f);
        } else {
            text(20.0f * u, ImVec2(x0 + 30.0f * u, 30.0f * u), IM_COL32(230, 220, 255, 240),
                 "Ejecta density, edge-on (slice through the spin axis)", 0.0f);
            text(16.0f * u, ImVec2(x0 + 30.0f * u, 56.0f * u), IM_COL32(190, 190, 210, 220),
                 "poles: fast, lanthanide-poor (blue)   plane: lanthanide-rich (red)", 0.0f);
            text(18.0f * u, ImVec2(x0 + w - 30.0f * u, 30.0f * u), IM_COL32(200, 200, 210, 230), "Time since the merger", 1.0f);
            text(34.0f * u, ImVec2(x0 + w - 30.0f * u, 54.0f * u), IM_COL32(255, 255, 255, 255), formatTime(t), 1.0f);
            // The kilonova's light curve: each component and their sum, log L against log t.
            const float bx = x0 + 30.0f * u, bw = std::min(380.0f * u, 0.55f * w), bh = 200.0f * u;
            const float by = H - bh - 40.0f * u;
            list.AddRectFilled(ImVec2(bx, by), ImVec2(bx + bw, by + bh), IM_COL32(0, 0, 0, 150));
            list.AddRect(ImVec2(bx, by), ImVec2(bx + bw, by + bh), IM_COL32(120, 120, 140, 200));
            const double lt0 = std::log10(0.03 * kDay), lt1 = std::log10(15.0 * kDay), lL0 = 39.5, lL1 = 42.5;
            auto px = [&](double tt, double L) {
                const float fx = static_cast<float>((std::log10(tt) - lt0) / (lt1 - lt0));
                const float fy = static_cast<float>((std::log10(std::max(L, 1.0)) - lL0) / (lL1 - lL0));
                return ImVec2(bx + bw * fx, by + bh * (1.0f - std::clamp(fy, 0.0f, 1.0f)));
            };
            const ImU32 cols[4] = {IM_COL32(90, 160, 255, 255), IM_COL32(230, 110, 160, 255), IM_COL32(255, 120, 50, 255),
                                   IM_COL32(255, 255, 255, 255)};
            for (int k = 0; k < 4; ++k) {
                std::vector<ImVec2> pts;
                for (int i = 0; i <= 80; ++i) {
                    const double tt = std::pow(10.0, lt0 + (lt1 - lt0) * i / 80.0);
                    double L = 0.0;
                    if (k < 3) L = m_kilonova[static_cast<size_t>(k)].luminosity(tt);
                    else
                        for (const Kilonova& kn : m_kilonova) L += kn.luminosity(tt);
                    pts.push_back(px(tt, L));
                }
                list.AddPolyline(pts.data(), static_cast<int>(pts.size()), cols[k], 0, (k == 3 ? 2.4f : 1.5f) * u);
            }
            if (t > 0.03 * kDay) {
                double L = 0.0;
                for (const Kilonova& kn : m_kilonova) L += kn.luminosity(t);
                const ImVec2 m = px(std::min(t, 15.0 * kDay), L);
                list.AddCircleFilled(m, 5.0f * u, IM_COL32(255, 255, 120, 255));
            }
            for (double d : {0.1, 1.0, 10.0}) {
                const ImVec2 a = px(d * kDay, 1e39);
                list.AddLine(ImVec2(a.x, by + bh), ImVec2(a.x, by + bh - 6.0f * u), IM_COL32(200, 200, 210, 220));
                std::snprintf(buf, sizeof(buf), d < 1.0 ? "%.1f d" : "%.0f d", d);
                text(14.0f * u, ImVec2(a.x, by + bh + 4.0f * u), IM_COL32(200, 200, 210, 220), buf, 0.5f);
            }
            text(15.0f * u, ImVec2(bx + 8.0f * u, by + 6.0f * u), IM_COL32(230, 230, 240, 240),
                 "Kilonova light curve (log L, 10^39.5 to 10^42.5 erg/s)", 0.0f);
        }
        // Scale bar.
        const double half = edgeOn ? 1.15 * ejectaFront(t) : m_viewHalf;
        const double km = niceKm(half * kKm * 0.5);
        const float px = static_cast<float>(km / kKm / (2.0 * half) * H);
        const float by = H - 50.0f * u, bx = x0 + w - 40.0f * u - px;
        list.AddLine(ImVec2(bx, by), ImVec2(bx + px, by), IM_COL32(255, 255, 255, 230), 3.0f * u);
        text(18.0f * u, ImVec2(bx + 0.5f * px, by - 30.0f * u), IM_COL32(255, 255, 255, 230), formatKm(km), 0.5f);
        // Phase.
        std::string phase = phaseName();
        text(20.0f * u, ImVec2(x0 + 0.5f * w, H - 60.0f * u), IM_COL32(220, 200, 255, 230), phase, 0.5f);
        // The color bar.
        const float cbx = x0 + w - 60.0f * u, cbTop = 160.0f * u, cbH = 260.0f * u;
        for (int k = 0; k < 32; ++k) {
            const float f0 = static_cast<float>(k) / 32.0f, f1 = static_cast<float>(k + 1) / 32.0f;
            // Inferno, roughly, for the legend.
            auto infer = [](float s) {
                const float r = std::clamp(-0.05f + 2.2f * s - 1.2f * s * s, 0.0f, 1.0f);
                const float g = std::clamp(s * s * 1.05f, 0.0f, 1.0f);
                const float b = std::clamp(0.1f + 1.6f * s - 2.6f * s * s + 1.2f * s * s * s, 0.0f, 1.0f);
                return IM_COL32(static_cast<int>(r * 255), static_cast<int>(g * 255), static_cast<int>(b * 255), 255);
            };
            list.AddRectFilled(ImVec2(cbx, cbTop + cbH * (1.0f - f1)), ImVec2(cbx + 14.0f * u, cbTop + cbH * (1.0f - f0)),
                               infer(0.5f * (f0 + f1)));
        }
        text(15.0f * u, ImVec2(cbx + 7.0f * u, cbTop + cbH + 8.0f * u), IM_COL32(200, 200, 210, 220),
             edgeOn ? "density" : "column", 0.5f);
        text(15.0f * u, ImVec2(cbx + 7.0f * u, cbTop + cbH + 26.0f * u), IM_COL32(200, 200, 210, 220),
             edgeOn ? "(log)" : "density", 0.5f);
        if (splitX > 0) list.AddLine(ImVec2(x0, 0.0f), ImVec2(x0, H), IM_COL32(60, 60, 70, 255), 2.0f);
    }
}

// The film's captions, as the reference sets them: bold italic, white with words picked out in
// yellow and red, at the top of the 9:16 frame; each new one comes in glowing with its letters
// spread, or sliding in, and the first leaves at the merger by glowing, spreading and smearing.
void MergerScene::fillFilmCaptions(ImDrawList& list, int width, int height) {
    ImFont* font = m_captionFont ? m_captionFont : ImGui::GetFont();
    const double t = m_filmTime;
    // The frame on the screen, and one of the reference's pixels.
    const float frameH = static_cast<float>(std::min(static_cast<double>(height), width / film::kAspect));
    const float u = frameH / static_cast<float>(film::kRefHeight);
    const float x0 = 0.5f * (static_cast<float>(width) - frameH * static_cast<float>(film::kAspect));
    const float y0 = 0.5f * (static_cast<float>(height) - frameH);
    const float frameW = frameH * static_cast<float>(film::kAspect);
    auto inkColour = [](film::Ink ink) {
        switch (ink) {
        case film::Ink::Yellow: return ImVec4(241 / 255.0f, 220 / 255.0f, 23 / 255.0f, 1.0f);
        case film::Ink::Red: return ImVec4(223 / 255.0f, 6 / 255.0f, 10 / 255.0f, 1.0f);
        default: return ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
        }
    };
    for (const film::Caption& c : film::captions()) {
        if (t < c.start || t >= c.end) continue;
        // The size that gives the caption its width in the reference (the widest line).
        // Segoe UI Black Italic at the reference's letter height (ImGui sizes a font by its
        // ascent plus descent, 1.33 of its em).
        const float nominal = 33.5f * u;
        float widest = 0.0f;
        for (const auto& line : c.lines) {
            std::string all;
            for (const film::Word& w : line) all += w.text;
            widest = std::max(widest, font->CalcTextSizeA(nominal, 1e9f, 0.0f, all.c_str()).x);
        }
        const float size = c.refWidth > 0.0 && widest > 0.0f
                               ? std::clamp(nominal * static_cast<float>(c.refWidth) * u / widest, 0.8f * nominal, 1.2f * nominal)
                               : nominal;
        ImFontBaked* baked = font->GetFontBaked(size);
        const float ascent = baked ? baked->Ascent : 0.8f * size;
        // How it comes in: e goes 1 -> 0 as it settles.
        const double since = t - c.start;
        float e = c.entry == film::Entry::None ? 0.0f : static_cast<float>(1.0 - smooth01(since / c.settle));
        float spread = 0.0f, glow = 0.0f, slide = 0.0f, smear = 0.0f, alpha = 1.0f, drop = 0.0f;
        if (c.entry == film::Entry::Glow) {
            spread = 0.35f * e;
            glow = e;
            alpha = 1.0f - 0.35f * e;
        } else if (c.entry == film::Entry::SlideFromLeft || c.entry == film::Entry::SlideFromRight) {
            slide = (c.entry == film::Entry::SlideFromLeft ? -1.0f : 1.0f) * 0.55f * frameW * e * e;
            glow = 0.8f * e;
            spread = 0.15f * e;
        }
        if (c.glitchOut) {
            // 6.13 - 6.55 s: a growing glow; 6.6 - 6.75: dimmer, sunk a little; 6.75 - 6.88:
            // letters spread wide; 6.88 - 6.95: smeared into streaks.
            glow = std::max(glow, static_cast<float>(smooth01((t - 6.1) / 0.2) * (1.0 - smooth01((t - 6.55) / 0.08))));
            if (t > 6.6) {
                alpha = 0.75f;
                drop = 6.0f * u;
            }
            if (t > 6.75) spread = 0.22f;
            if (t > 6.88) smear = 1.0f;
        }
        // Its lines, laid out letter by letter so they can spread.
        const bool twoLines = c.lines.size() > 1;
        for (size_t li = 0; li < c.lines.size(); ++li) {
            const auto& line = c.lines[li];
            std::string all;
            for (const film::Word& w : line) all += w.text;
            const float plain = font->CalcTextSizeA(size, 1e9f, 0.0f, all.c_str()).x;
            int letters = 0;
            for (const film::Word& w : line) letters += static_cast<int>(w.text.size());
            const float gap = spread * size * 0.5f;
            const float lineW = plain + gap * static_cast<float>(std::max(letters - 1, 0));
            // Baselines: one line at 59 reference pixels; the title's at 43 and 75.
            const float baseline = y0 + u * (twoLines ? (li == 0 ? 43.0f : 75.0f) : 59.0f) + drop;
            float x = c.leftAligned ? x0 + 13.0f * u : x0 + 0.5f * frameW - 0.5f * lineW;
            x += slide;
            const float top = baseline - ascent;
            for (const film::Word& w : line) {
                const ImVec4 ink = inkColour(w.ink);
                for (size_t ci = 0; ci < w.text.size(); ++ci) {
                    const char* ch = w.text.c_str() + ci;
                    const float adv = font->CalcTextSizeA(size, 1e9f, 0.0f, ch, ch + 1).x;
                    const ImVec2 at(x, top);
                    // A soft dark shadow, so it reads over the bright shots.
                    for (int k = 0; k < 4; ++k) {
                        const float dx = (k & 1 ? 1.0f : -1.0f) * 1.2f * u, dy = (k & 2 ? 1.0f : -1.0f) * 1.2f * u;
                        list.AddText(font, size, ImVec2(at.x + dx + 1.0f * u, at.y + dy + 1.5f * u),
                                     ImGui::GetColorU32(ImVec4(0.0f, 0.0f, 0.0f, 0.22f * alpha)), ch, ch + 1);
                    }
                    // The glow: blurred, brighter copies around the letter.
                    if (glow > 0.01f) {
                        const ImVec4 g(std::min(ink.x + 0.3f, 1.0f), std::min(ink.y + 0.3f, 1.0f), std::min(ink.z + 0.3f, 1.0f),
                                       0.16f * glow);
                        for (int k = 0; k < 8; ++k) {
                            const float a = static_cast<float>(k) * 0.785398f;
                            const float r = 3.5f * u * (1.0f + glow);
                            list.AddText(font, size, ImVec2(at.x + r * std::cos(a), at.y + r * std::sin(a)), ImGui::GetColorU32(g),
                                         ch, ch + 1);
                        }
                    }
                    if (smear > 0.0f) {
                        for (int k = -8; k <= 8; ++k)
                            list.AddText(font, size, ImVec2(at.x, at.y + 2.6f * u * static_cast<float>(k)),
                                         ImGui::GetColorU32(ImVec4(ink.x, ink.y, ink.z, 0.09f)), ch, ch + 1);
                    } else {
                        list.AddText(font, size, at, ImGui::GetColorU32(ImVec4(ink.x, ink.y, ink.z, alpha)), ch, ch + 1);
                    }
                    // Sliding in, it trails a motion blur.
                    if (std::abs(slide) > 1.0f) {
                        for (int k = 1; k <= 4; ++k)
                            list.AddText(font, size, ImVec2(at.x - slide * 0.06f * static_cast<float>(k), at.y),
                                         ImGui::GetColorU32(ImVec4(ink.x, ink.y, ink.z, 0.18f)), ch, ch + 1);
                    }
                    x += adv + gap;
                }
            }
        }
    }
}

std::string MergerScene::phaseName() const {
    if (view == Film) return std::string("Film: ") + film::shotName(m_filmTime);
    const double t = sinceMerger();
    if (!m_sim.merged()) return m_sim.phase() == MergerSim::Phase::Inspiral ? "Inspiral" : "Inspiral (hydrodynamics)";
    if (!m_sim.hasHole()) return "Hypermassive neutron star";
    if (m_sim.phase() != MergerSim::Phase::Coasting) return "Black hole and torus";
    if (jetOn && t > jet.launch && t < burstTime() + 6.0) return "Jet and gamma-ray burst";
    if (t < 3600.0) return "Ejecta coasting";
    return "Kilonova";
}

std::string MergerScene::describe() const {
    char buf[512];
    if (view == Film) {
        std::snprintf(buf, sizeof(buf), "merger film: t = %.2f s of %.1f, %s", m_filmTime, film::kLength,
                      film::shotName(m_filmTime));
        return buf;
    }
    std::snprintf(buf, sizeof(buf),
                  "merger: %s, t = %s (sim %.3f ms, %d steps, last dt %.3f), separation %.1f km, rho_max %.3g g/cc, "
                  "T_max %.3g K, unbound %.4f Msun%s",
                  phaseName().c_str(), formatTime(sinceMerger()).c_str(), m_sim.time() * kSecond * 1e3, m_sim.stepsTaken(),
                  m_sim.lastStep(), m_sim.separation() * kKm, m_sim.maxDensity(), m_sim.maxTemperature(), m_sim.unboundMass(),
                  m_sim.hasHole() ? (", hole " + std::to_string(m_sim.holeMass()) + " Msun").c_str() : "");
    char cam[160];
    std::snprintf(cam, sizeof(cam), "; camera at %s (%s), blue ejecta out to %s", formatKm(m_lastR * kKm).c_str(),
                  autoCamera ? "following" : "held", formatKm(ejectaFront(sinceMerger()) * kKm).c_str());
    return std::string(buf) + cam;
}

// --- Panel ------------------------------------------------------------------------------------------------------

int MergerScene::drawControls(Camera& camera) {
    int request = 0;
    if (!m_active) {
        ImGui::TextWrapped("Two neutron stars spiral together, merge, collapse into a black hole, fire a gamma-ray "
                           "burst and glow as a kilonova, as GW170817 did in 2017: smoothed-particle hydrodynamics "
                           "on the GPU, then the kilonova's light curves.");
        float m1 = static_cast<float>(params.m1), m2 = static_cast<float>(params.m2);
        float radius = static_cast<float>(params.radiusKm), collapse = static_cast<float>(params.collapseMs);
        if (ImGui::SliderFloat("Star 1 (Suns)", &m1, 1.1f, 1.6f, "%.2f")) params.m1 = m1;
        if (ImGui::SliderFloat("Star 2 (Suns)", &m2, 1.1f, 1.6f, "%.2f")) params.m2 = m2;
        if (ImGui::SliderFloat("Radius (km)", &radius, 10.0f, 14.0f, "%.1f")) params.radiusKm = radius;
        if (ImGui::SliderFloat("Collapse after (ms)", &collapse, 3.0f, 40.0f, "%.0f")) params.collapseMs = collapse;
        ImGui::SetItemTooltip("How long the hypermassive neutron star lasts before it collapses into a black hole.");
        static const int kResolutions[] = {4000, 8000, 12000, 20000};
        static const char* kResolutionNames[] = {"4,000 per star (fast)", "8,000 per star", "12,000 per star",
                                                 "20,000 per star (slow)"};
        int res = 2;
        for (int k = 0; k < 4; ++k)
            if (params.perStar == kResolutions[k]) res = k;
        if (ImGui::Combo("Particles", &res, kResolutionNames, 4)) params.perStar = kResolutions[res];
        if (ImGui::Button("Start the merger", ImVec2(-1.0f, 0.0f))) request = 1;
        return request;
    }
    ImGui::TextUnformatted("View");
    ImGui::RadioButton("Cinematic", &view, Cinematic);
    ImGui::SameLine();
    ImGui::RadioButton("Simulation", &view, Simulation);
    ImGui::SameLine();
    ImGui::RadioButton("Side by side", &view, SideBySide);
    ImGui::SameLine();
    if (ImGui::RadioButton("Film", &view, Film)) setFilmTime(m_filmTime);
    ImGui::SetItemTooltip("The reference film's shots (GW170817 as NASA animated it), replayed with its timing, "
                          "framing, colours and captions. It is framed 9:16: record it at 1080 x 1920.");
    ImGui::SliderFloat("Speed", &speed, 0.1f, 5.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);
    if (view == Film) {
        float ft = static_cast<float>(m_filmTime);
        if (ImGui::SliderFloat("Film time", &ft, 0.0f, static_cast<float>(film::kLength) - 0.01f, "%.2f s"))
            setFilmTime(ft);
        ImGui::Checkbox("Captions", &captions);
        ImGui::SameLine();
        ImGui::Checkbox("Sound", &sound);
        if (sound) ImGui::SliderFloat("Volume", &volume, 0.0f, 1.0f, "%.2f");
        ImGui::Separator();
        ImGui::TextWrapped("%s", film::shotName(m_filmTime));
        if (ImGui::Button("Restart")) setFilmTime(0.0);
        ImGui::SameLine();
        if (ImGui::Button("Leave the merger")) request = 2;
        return request;
    }
    ImGui::Checkbox("Captions", &captions);
    ImGui::SameLine();
    ImGui::Checkbox("Gravitational waves", &waves);
    ImGui::Checkbox("Debris cloud", &cloudOn);
    ImGui::SameLine();
    ImGui::Checkbox("Kilonova", &kilonovaOn);
    ImGui::SameLine();
    ImGui::Checkbox("Jet and burst", &jetOn);
    ImGui::Checkbox("Film colors##merger", &filmColor);
    ImGui::SameLine();
    if (ImGui::Checkbox("Follow the action", &autoCamera) && autoCamera) {
        m_autoAngles = true;
        m_lastR = camera.r;
    }
    ImGui::Checkbox("Sound", &sound);
    if (sound) {
        ImGui::SliderFloat("Volume", &volume, 0.0f, 1.0f, "%.2f");
        ImGui::SliderFloat("Pitch", &pitch, 0.25f, 2.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);
        ImGui::SetItemTooltip("1x plays the gravitational wave at its true frequency: it is audible.");
    }
    ImGui::Separator();
    const double t = sinceMerger();
    ImGui::Text("%s, t = %s", phaseName().c_str(), formatTime(t).c_str());
    ImGui::Text("Simulated %.3f ms in %d steps (%.1f ms GPU)", m_sim.time() * kSecond * 1e3, m_sim.stepsTaken(),
                m_sim.gpuMilliseconds());
    if (!m_sim.merged()) {
        // Twice the orbital frequency (Kepler's, at the stars' separation).
        const double a = m_sim.separation();
        ImGui::Text("Separation %.1f km, waves at %.0f Hz", a * kKm,
                    std::sqrt(m_sim.totalMass() / (a * a * a)) / kPi / kSecond);
    }
    ImGui::Text("Densest %.2g g/cm3, hottest %.2g K", m_sim.maxDensity(), m_sim.maxTemperature());
    ImGui::Text("Unbound %.4f Suns", m_sim.unboundMass());
    if (m_sim.hasHole())
        ImGui::Text("Black hole %.3f Suns, spin %.2f, torus %.3f Suns", m_sim.holeMass(), m_sim.hole().spin,
                    m_sim.hole().diskMass);
    if (t > 0.01) {
        double L = 0.0;
        for (const Kilonova& k : m_kilonova) L += k.luminosity(t);
        ImGui::Text("Kilonova %.2g erg/s", L);
    }
    if (ImGui::Button("Restart")) request = 1;
    ImGui::SameLine();
    if (ImGui::Button("Leave the merger")) request = 2;
    return request;
}
