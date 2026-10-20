#pragma once

// The scene's sound. Space carries none, so this is a film's convention, but every sound
// is carried by light: its pitch and its rate are multiplied by g, the same frequency ratio
// the ray tracer colours the light with (gravitational redshift and Doppler shift).
//   bodies     a soft tone each, pitched by kind and temperature, louder close by, panned
//              by where the body is in the view;
//   pulsars    a click each time a beam sweeps past (twice a turn), at the shifted rate,
//              louder by g^3 (beaming);
//   rumble     a low roar that grows with the tide at the camera, ~ 2M / r^3;
//   waves      heavy bodies send out gravitational waves at twice their orbital frequency,
//              raised in pitch into hearing (the true frequency is shown).
// Played in stereo (waveOut); videos get a mono mix.

#include "audio.h"

#include <cstdint>
#include <map>
#include <vector>

class SceneSound {
public:
    static constexpr int kSampleRate = 44100;
    enum Kind { Tone, Clicks, Wave };
    struct Voice {
        int id = 0;               // keeps each voice's phase from frame to frame
        Kind kind = Tone;
        double frequency = 220.0; // Hz as heard (for clicks: of the ping)
        double gain = 0.0;        // 0 to 1
        double pan = 0.0;         // -1 left to 1 right
        double rate = 0.0;        // clicks per second as heard
    };
    // Makes wallDt seconds of sound (live: plays it, keeping ~0.12 s queued).
    void update(const std::vector<Voice>& voices, double rumble, double wallDt, bool live, bool enabled, double volume);
    // The sound made by the last update, mono, for videos.
    const std::vector<int16_t>& lastMono() const { return m_mono; }
    void close() { m_audio.close(); }

private:
    struct State {
        double phase = 0.0, frequency = 0.0, gain = 0.0, pan = 0.0;
        double clickPhase = 0.0, clickAge = 1e9;
        bool started = false;
    };
    std::map<int, State> m_states;
    AudioStream m_audio;
    double m_carry = 0.0;
    double m_rumble = 0.0, m_low1 = 0.0, m_low2 = 0.0;
    uint32_t m_noise = 12345u;
    std::vector<int16_t> m_mono;
};
