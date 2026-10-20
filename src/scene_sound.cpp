#include "scene_sound.h"

#include <algorithm>
#include <cmath>

namespace {

constexpr double kTwoPi = 6.283185307179586;
constexpr double kQueued = 0.12;        // seconds of sound kept waiting to play
constexpr double kClickDecay = 0.004;   // seconds: a pulsar's tick
constexpr double kClickPing = 2600.0;   // Hz: what it rings at (times g)
constexpr double kRumbleCutoff = 0.01;  // the rumble's low-pass, per sample (~70 Hz)

}  // namespace

void SceneSound::update(const std::vector<Voice>& voices, double rumble, double wallDt, bool live, bool enabled,
                        double volume) {
    m_mono.clear();
    if (live && enabled && !m_audio.isOpen()) m_audio.open(kSampleRate, 2);
    if ((!live || !enabled) && m_audio.isOpen()) m_audio.close();
    if (!enabled) {
        m_states.clear();
        m_carry = 0.0;
        return;
    }
    double want = wallDt * kSampleRate + m_carry;
    if (live) want += (kQueued - m_audio.queued()) * kSampleRate * 0.5;
    const int n = std::max(0, static_cast<int>(want));
    m_carry = live ? 0.0 : want - n;
    if (n == 0) return;

    // Voices that went silent fade out over this piece, then are forgotten.
    for (auto& [id, st] : m_states) st.started = false;
    struct Ramp {
        State* st;
        const Voice* v;
        double f0, f1, g0, g1, p0, p1;
    };
    std::vector<Ramp> ramps;
    for (const Voice& v : voices) {
        State& st = m_states[v.id];
        if (st.frequency <= 0.0) {
            st.frequency = v.frequency;
            st.pan = v.pan;
        }
        st.started = true;
        ramps.push_back({&st, &v, st.frequency, v.frequency, st.gain, v.gain, st.pan, v.pan});
    }
    static const Voice silent{};
    for (auto& [id, st] : m_states) {
        if (!st.started) ramps.push_back({&st, &silent, st.frequency, st.frequency, st.gain, 0.0, st.pan, st.pan});
    }

    std::vector<int16_t> stereo(static_cast<size_t>(2 * n));
    m_mono.resize(static_cast<size_t>(n));
    const double rumble0 = m_rumble;
    for (int k = 0; k < n; ++k) {
        const double u = (k + 1.0) / n;
        double left = 0.0, right = 0.0;
        for (const Ramp& r : ramps) {
            State& st = *r.st;
            const double f = r.f0 + (r.f1 - r.f0) * u;
            const double gain = r.g0 + (r.g1 - r.g0) * u;
            const double pan = r.p0 + (r.p1 - r.p0) * u;
            double s = 0.0;
            if (r.v->kind == Clicks) {
                // A tick each time a beam sweeps past.
                st.clickPhase += r.v->rate / kSampleRate;
                if (st.clickPhase >= 1.0) {
                    st.clickPhase -= std::floor(st.clickPhase);
                    st.clickAge = 0.0;
                }
                if (st.clickAge < 10.0 * kClickDecay) {
                    s = gain * std::exp(-st.clickAge / kClickDecay) * std::sin(kTwoPi * f * st.clickAge);
                }
                st.clickAge += 1.0 / kSampleRate;
            } else {
                st.phase = std::fmod(st.phase + kTwoPi * f / kSampleRate, kTwoPi);
                if (r.v->kind == Wave) {
                    s = gain * std::sin(st.phase);
                } else {
                    // A soft tone: mostly the fundamental, a little of the next two harmonics.
                    s = gain * (0.78 * std::sin(st.phase) + 0.17 * std::sin(2.0 * st.phase) + 0.05 * std::sin(3.0 * st.phase));
                }
            }
            // Equal-power panning.
            const double angle = 0.25 * kTwoPi * 0.5 * (std::clamp(pan, -1.0, 1.0) + 1.0);
            left += s * std::cos(angle);
            right += s * std::sin(angle);
        }
        // The rumble: white noise through two low-pass stages.
        const double level = rumble0 + (rumble - rumble0) * u;
        if (level > 1e-4) {
            m_noise = m_noise * 1664525u + 1013904223u;
            const double white = (m_noise >> 9) / 4194304.0 - 1.0;
            m_low1 += (white - m_low1) * kRumbleCutoff;
            m_low2 += (m_low1 - m_low2) * kRumbleCutoff;
            const double roar = 9.0 * level * m_low2;
            left += roar;
            right += roar;
        }
        // Gentle limiting, so many voices together never clip harshly.
        left = std::tanh(volume * left);
        right = std::tanh(volume * right);
        stereo[static_cast<size_t>(2 * k)] = static_cast<int16_t>(left * 30000.0);
        stereo[static_cast<size_t>(2 * k + 1)] = static_cast<int16_t>(right * 30000.0);
        m_mono[static_cast<size_t>(k)] = static_cast<int16_t>(0.5 * (left + right) * 30000.0);
    }
    m_rumble = rumble;
    for (const Ramp& r : ramps) {
        r.st->frequency = r.f1;
        r.st->gain = r.g1;
        r.st->pan = r.p1;
    }
    for (auto it = m_states.begin(); it != m_states.end();) {
        if (!it->second.started && it->second.gain <= 0.0) it = m_states.erase(it);
        else ++it;
    }
    if (live) m_audio.submit(stereo);
}
