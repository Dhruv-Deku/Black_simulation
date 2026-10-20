#pragma once

#include <cstdint>
#include <memory>
#include <vector>

// Plays a stream of 16-bit samples (mono, or interleaved stereo) through Windows' waveOut (part of Windows, so no
// extra software is needed). The app synthesizes a little more each frame; queued() says
// how much is waiting so it can keep about a tenth of a second ahead.
class AudioStream {
public:
    AudioStream();
    ~AudioStream();
    AudioStream(const AudioStream&) = delete;
    AudioStream& operator=(const AudioStream&) = delete;

    bool open(int sampleRate, int channels = 1);
    void close();
    bool isOpen() const;
    int sampleRate() const { return m_rate; }
    void submit(const std::vector<int16_t>& samples);
    double queued() const;  // seconds of sound waiting to be played

private:
    struct Impl;
    std::unique_ptr<Impl> m;
    int m_rate = 44100;
    int m_channels = 1;
};
