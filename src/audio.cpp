#include "audio.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <mmsystem.h>
#endif

#include <algorithm>
#include <cstring>

#ifdef _WIN32
struct AudioStream::Impl {
    static constexpr int kBuffers = 16;
    HWAVEOUT device = nullptr;
    WAVEHDR headers[kBuffers] = {};
    std::vector<int16_t> data[kBuffers];
    long long submitted = 0;  // samples handed to the device
};

AudioStream::AudioStream() : m(std::make_unique<Impl>()) {}
AudioStream::~AudioStream() {
    close();
}

bool AudioStream::open(int sampleRate, int channels) {
    close();
    m_rate = sampleRate;
    m_channels = channels;
    WAVEFORMATEX format{};
    format.wFormatTag = WAVE_FORMAT_PCM;
    format.nChannels = static_cast<WORD>(channels);
    format.nSamplesPerSec = static_cast<DWORD>(sampleRate);
    format.wBitsPerSample = 16;
    format.nBlockAlign = static_cast<WORD>(2 * channels);
    format.nAvgBytesPerSec = format.nSamplesPerSec * format.nBlockAlign;
    if (waveOutOpen(&m->device, WAVE_MAPPER, &format, 0, 0, CALLBACK_NULL) != MMSYSERR_NOERROR) {
        m->device = nullptr;
        return false;
    }
    m->submitted = 0;
    return true;
}

void AudioStream::close() {
    if (!m || !m->device) return;
    waveOutReset(m->device);
    for (WAVEHDR& h : m->headers) {
        if (h.dwFlags & WHDR_PREPARED) waveOutUnprepareHeader(m->device, &h, sizeof(WAVEHDR));
        h = WAVEHDR{};
    }
    waveOutClose(m->device);
    m->device = nullptr;
}

bool AudioStream::isOpen() const {
    return m && m->device;
}

void AudioStream::submit(const std::vector<int16_t>& samples) {
    if (!isOpen() || samples.empty()) return;
    for (int k = 0; k < Impl::kBuffers; ++k) {
        WAVEHDR& h = m->headers[k];
        const bool busy = (h.dwFlags & WHDR_PREPARED) && !(h.dwFlags & WHDR_DONE);
        if (busy) continue;
        if (h.dwFlags & WHDR_PREPARED) waveOutUnprepareHeader(m->device, &h, sizeof(WAVEHDR));
        m->data[k] = samples;
        h = WAVEHDR{};
        h.lpData = reinterpret_cast<LPSTR>(m->data[k].data());
        h.dwBufferLength = static_cast<DWORD>(m->data[k].size() * sizeof(int16_t));
        waveOutPrepareHeader(m->device, &h, sizeof(WAVEHDR));
        waveOutWrite(m->device, &h, sizeof(WAVEHDR));
        m->submitted += static_cast<long long>(samples.size()) / m_channels;  // sample frames
        return;
    }
    // Every buffer is still playing: drop this piece rather than wait.
}

double AudioStream::queued() const {
    if (!isOpen()) return 0.0;
    MMTIME t{};
    t.wType = TIME_SAMPLES;
    waveOutGetPosition(m->device, &t, sizeof(MMTIME));
    const long long played = t.wType == TIME_SAMPLES ? static_cast<long long>(t.u.sample) : m->submitted;
    return std::max(0.0, static_cast<double>(m->submitted - played) / m_rate);
}
#else
struct AudioStream::Impl {};
AudioStream::AudioStream() : m(std::make_unique<Impl>()) {}
AudioStream::~AudioStream() = default;
bool AudioStream::open(int sampleRate, int channels) {
    m_rate = sampleRate;
    m_channels = channels;
    return false;
}
void AudioStream::close() {}
bool AudioStream::isOpen() const {
    return false;
}
void AudioStream::submit(const std::vector<int16_t>&) {}
double AudioStream::queued() const {
    return 0.0;
}
#endif
