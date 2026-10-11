#pragma once

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

    // width and height must be even.
    bool open(const std::string& path, int width, int height, int fps, std::string& error);
    // One frame of top-down BGRA bytes (width * height * 4).
    bool write(const unsigned char* bgra, std::string& error);
    bool close(std::string& error);
    bool isOpen() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m;
};
