#include "recorder.h"

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <codecapi.h>
#include <mfapi.h>
#include <mferror.h>
#include <mfidl.h>
#include <mfreadwrite.h>
#include <wrl/client.h>

#include <algorithm>
#include <cstdio>
#include <cstring>

using Microsoft::WRL::ComPtr;

struct VideoRecorder::Impl {
    ComPtr<IMFSinkWriter> writer;
    DWORD stream = 0;
    DWORD audioStream = 0;
    bool audio = false;
    int audioRate = 0;
    long long audioSamples = 0;
    int width = 0, height = 0, fps = 30;
    long long frame = 0;
    bool comInitialized = false;
    bool mfStarted = false;
};

namespace {

std::string failure(const char* what, HRESULT hr) {
    char buf[160];
    std::snprintf(buf, sizeof(buf), "%s failed (HRESULT 0x%08lX)", what, static_cast<unsigned long>(hr));
    return buf;
}

std::wstring widen(const std::string& s) {
    int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, nullptr, 0);
    std::wstring w(static_cast<size_t>(std::max(n, 1)), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, w.data(), n);
    w.resize(wcslen(w.c_str()));
    return w;
}

HRESULT setVideoFormat(IMFMediaType* type, int width, int height, int fps) {
    HRESULT hr = type->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
    if (SUCCEEDED(hr)) hr = type->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive);
    if (SUCCEEDED(hr)) hr = MFSetAttributeSize(type, MF_MT_FRAME_SIZE, width, height);
    if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(type, MF_MT_FRAME_RATE, fps, 1);
    if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(type, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);
    return hr;
}

}  // namespace

VideoRecorder::VideoRecorder() : m(std::make_unique<Impl>()) {}

VideoRecorder::~VideoRecorder() {
    std::string ignored;
    close(ignored);
}

bool VideoRecorder::isOpen() const {
    return m->writer != nullptr;
}

bool VideoRecorder::open(const std::string& path, int width, int height, int fps, std::string& error, int audioRate) {
    if (isOpen()) close(error);
    if (width % 2 || height % 2 || width <= 0 || height <= 0) {
        error = "video width and height must be positive and even";
        return false;
    }
    // S_FALSE means COM was already initialized on this thread; still balance it.
    HRESULT hr = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
    m->comInitialized = SUCCEEDED(hr);
    hr = MFStartup(MF_VERSION);
    if (FAILED(hr)) {
        error = failure("Starting Media Foundation", hr);
        return false;
    }
    m->mfStarted = true;

    ComPtr<IMFAttributes> attributes;
    MFCreateAttributes(&attributes, 2);
    attributes->SetUINT32(MF_READWRITE_ENABLE_HARDWARE_TRANSFORMS, TRUE);  // GPU encoder when available
    attributes->SetUINT32(MF_SINK_WRITER_DISABLE_THROTTLING, TRUE);
    hr = MFCreateSinkWriterFromURL(widen(path).c_str(), nullptr, attributes.Get(), &m->writer);
    if (FAILED(hr)) {
        error = failure("Creating the video file", hr);
        close(error);
        return false;
    }

    // Output: H.264 High profile at about 0.12 bits per pixel per frame.
    const UINT32 bitrate = static_cast<UINT32>(std::clamp(0.12 * width * height * fps, 4e6, 100e6));
    ComPtr<IMFMediaType> out;
    MFCreateMediaType(&out);
    hr = setVideoFormat(out.Get(), width, height, fps);
    if (SUCCEEDED(hr)) hr = out->SetGUID(MF_MT_SUBTYPE, MFVideoFormat_H264);
    if (SUCCEEDED(hr)) hr = out->SetUINT32(MF_MT_AVG_BITRATE, bitrate);
    if (SUCCEEDED(hr)) hr = out->SetUINT32(MF_MT_MPEG2_PROFILE, eAVEncH264VProfile_High);
    if (SUCCEEDED(hr)) hr = m->writer->AddStream(out.Get(), &m->stream);
    if (FAILED(hr)) {
        error = failure("Setting up the H.264 encoder", hr);
        m->writer.Reset();
        close(error);
        return false;
    }

    // Input: 32-bit BGRX rows, top-down (a positive stride); the sink writer
    // converts to the encoder's YUV format itself.
    ComPtr<IMFMediaType> in;
    MFCreateMediaType(&in);
    hr = setVideoFormat(in.Get(), width, height, fps);
    if (SUCCEEDED(hr)) hr = in->SetGUID(MF_MT_SUBTYPE, MFVideoFormat_RGB32);
    if (SUCCEEDED(hr)) hr = in->SetUINT32(MF_MT_DEFAULT_STRIDE, static_cast<UINT32>(width * 4));
    if (SUCCEEDED(hr)) hr = m->writer->SetInputMediaType(m->stream, in.Get(), nullptr);

    // Sound (optional): mono 16-bit PCM in, AAC at 128 kbit/s out. Without an AAC encoder
    // the video is still written, silent.
    m->audio = false;
    m->audioSamples = 0;
    if (SUCCEEDED(hr) && audioRate > 0) {
        ComPtr<IMFMediaType> aout, ain;
        MFCreateMediaType(&aout);
        HRESULT ha = aout->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Audio);
        if (SUCCEEDED(ha)) ha = aout->SetGUID(MF_MT_SUBTYPE, MFAudioFormat_AAC);
        if (SUCCEEDED(ha)) ha = aout->SetUINT32(MF_MT_AUDIO_BITS_PER_SAMPLE, 16);
        if (SUCCEEDED(ha)) ha = aout->SetUINT32(MF_MT_AUDIO_SAMPLES_PER_SECOND, static_cast<UINT32>(audioRate));
        if (SUCCEEDED(ha)) ha = aout->SetUINT32(MF_MT_AUDIO_NUM_CHANNELS, 1);
        if (SUCCEEDED(ha)) ha = aout->SetUINT32(MF_MT_AUDIO_AVG_BYTES_PER_SECOND, 16000);
        DWORD index = 0;
        if (SUCCEEDED(ha)) ha = m->writer->AddStream(aout.Get(), &index);
        MFCreateMediaType(&ain);
        if (SUCCEEDED(ha)) ha = ain->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Audio);
        if (SUCCEEDED(ha)) ha = ain->SetGUID(MF_MT_SUBTYPE, MFAudioFormat_PCM);
        if (SUCCEEDED(ha)) ha = ain->SetUINT32(MF_MT_AUDIO_BITS_PER_SAMPLE, 16);
        if (SUCCEEDED(ha)) ha = ain->SetUINT32(MF_MT_AUDIO_SAMPLES_PER_SECOND, static_cast<UINT32>(audioRate));
        if (SUCCEEDED(ha)) ha = ain->SetUINT32(MF_MT_AUDIO_NUM_CHANNELS, 1);
        if (SUCCEEDED(ha)) ha = ain->SetUINT32(MF_MT_AUDIO_BLOCK_ALIGNMENT, 2);
        if (SUCCEEDED(ha)) ha = ain->SetUINT32(MF_MT_AUDIO_AVG_BYTES_PER_SECOND, static_cast<UINT32>(audioRate * 2));
        if (SUCCEEDED(ha)) ha = m->writer->SetInputMediaType(index, ain.Get(), nullptr);
        if (SUCCEEDED(ha)) {
            m->audio = true;
            m->audioStream = index;
            m->audioRate = audioRate;
        } else {
            std::fprintf(stderr, "No sound in the video: %s\n", failure("Setting up the AAC encoder", ha).c_str());
        }
    }
    if (SUCCEEDED(hr)) hr = m->writer->BeginWriting();
    if (FAILED(hr)) {
        error = failure("Configuring the video input", hr);
        m->writer.Reset();
        close(error);
        return false;
    }
    m->width = width;
    m->height = height;
    m->fps = fps;
    m->frame = 0;
    return true;
}

bool VideoRecorder::write(const unsigned char* bgra, std::string& error) {
    if (!isOpen()) {
        error = "no video is being recorded";
        return false;
    }
    const DWORD size = static_cast<DWORD>(m->width) * m->height * 4;
    ComPtr<IMFMediaBuffer> buffer;
    HRESULT hr = MFCreateMemoryBuffer(size, &buffer);
    BYTE* dst = nullptr;
    if (SUCCEEDED(hr)) hr = buffer->Lock(&dst, nullptr, nullptr);
    if (SUCCEEDED(hr)) {
        std::memcpy(dst, bgra, size);
        buffer->Unlock();
        hr = buffer->SetCurrentLength(size);
    }
    ComPtr<IMFSample> sample;
    if (SUCCEEDED(hr)) hr = MFCreateSample(&sample);
    if (SUCCEEDED(hr)) hr = sample->AddBuffer(buffer.Get());
    const LONGLONG ticks = 10'000'000LL;  // Media Foundation time is in 100 ns units
    if (SUCCEEDED(hr)) hr = sample->SetSampleTime(m->frame * ticks / m->fps);
    if (SUCCEEDED(hr)) hr = sample->SetSampleDuration(ticks / m->fps);
    if (SUCCEEDED(hr)) hr = m->writer->WriteSample(m->stream, sample.Get());
    if (FAILED(hr)) {
        error = failure("Encoding a frame", hr);
        return false;
    }
    ++m->frame;
    return true;
}

bool VideoRecorder::writeAudio(const int16_t* samples, size_t count, std::string& error) {
    if (!isOpen() || !m->audio || count == 0) return true;
    const DWORD size = static_cast<DWORD>(count * sizeof(int16_t));
    ComPtr<IMFMediaBuffer> buffer;
    HRESULT hr = MFCreateMemoryBuffer(size, &buffer);
    BYTE* dst = nullptr;
    if (SUCCEEDED(hr)) hr = buffer->Lock(&dst, nullptr, nullptr);
    if (SUCCEEDED(hr)) {
        std::memcpy(dst, samples, size);
        buffer->Unlock();
        hr = buffer->SetCurrentLength(size);
    }
    ComPtr<IMFSample> sample;
    if (SUCCEEDED(hr)) hr = MFCreateSample(&sample);
    if (SUCCEEDED(hr)) hr = sample->AddBuffer(buffer.Get());
    const LONGLONG ticks = 10'000'000LL;
    if (SUCCEEDED(hr)) hr = sample->SetSampleTime(m->audioSamples * ticks / m->audioRate);
    if (SUCCEEDED(hr)) hr = sample->SetSampleDuration(static_cast<LONGLONG>(count) * ticks / m->audioRate);
    if (SUCCEEDED(hr)) hr = m->writer->WriteSample(m->audioStream, sample.Get());
    if (FAILED(hr)) {
        error = failure("Encoding sound", hr);
        return false;
    }
    m->audioSamples += static_cast<long long>(count);
    return true;
}

bool VideoRecorder::close(std::string& error) {
    bool ok = true;
    if (m->writer) {
        HRESULT hr = m->writer->Finalize();
        if (FAILED(hr)) {
            error = failure("Finishing the video file", hr);
            ok = false;
        }
        m->writer.Reset();
    }
    if (m->mfStarted) {
        MFShutdown();
        m->mfStarted = false;
    }
    if (m->comInitialized) {
        CoUninitialize();
        m->comInitialized = false;
    }
    return ok;
}

#else  // not Windows

struct VideoRecorder::Impl {};
VideoRecorder::VideoRecorder() : m(std::make_unique<Impl>()) {}
VideoRecorder::~VideoRecorder() = default;
bool VideoRecorder::isOpen() const { return false; }
bool VideoRecorder::open(const std::string&, int, int, int, std::string& error, int) {
    error = "video recording needs Windows Media Foundation";
    return false;
}
bool VideoRecorder::write(const unsigned char*, std::string& error) {
    error = "no video is being recorded";
    return false;
}
bool VideoRecorder::writeAudio(const int16_t*, size_t, std::string&) { return true; }
bool VideoRecorder::close(std::string&) { return true; }

#endif
