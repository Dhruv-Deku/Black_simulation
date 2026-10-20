#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

// Writes frames to an H.264 MP4 file using Windows Media Foundation, which ships
// with Windows, so no extra software is needed.
class VideoRecorder {
public:
    VideoRecorder();
    ~VideoRecorder();
    VideoRecorder(const VideoRecorder&) = delete;
    VideoRecorder& operator=(const VideoRecorder&) = delete;

    // width and height must be even. audioRate > 0 adds a mono AAC sound track fed by
    // writeAudio (16-bit samples at that rate).
    bool open(const std::string& path, int width, int height, int fps, std::string& error, int audioRate = 0);
    // One frame of top-down BGRA bytes (width * height * 4).
    bool write(const unsigned char* bgra, std::string& error);
    bool writeAudio(const int16_t* samples, size_t count, std::string& error);
    bool close(std::string& error);
    bool isOpen() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m;
};
