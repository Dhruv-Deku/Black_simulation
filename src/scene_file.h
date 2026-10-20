#pragma once

// Saving and loading scenes.
//
// A scene is two files: NAME.bhscene, readable JSON with everything that decides what
// happens next (the central object, every setting, the camera and its geodesic, every
// body with its constants of motion, state, shape and past, the bodies that are gone
// but whose light is still on its way, the tidal disruptions so far, the events, the
// time and the random number generators' states), and NAME.bhscene.bin beside it, the
// bulky binary parts: the light-delay history of every body and the debris gas (the
// app adds that; see App::saveScene).
//
// Doubles and floats are written in their shortest exact form, so loading a scene
// restores it bit for bit: saving it again gives the same files, and it runs on exactly
// as the original would have.

#include "bodies.h"
#include "camera.h"
#include "json.h"
#include "settings.h"

#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace scene {

constexpr int kVersion = 1;
constexpr const char* kExtension = ".bhscene";

// The binary sidecar: named chunks of bytes.
class Blob {
public:
    // Adds a chunk (tags are up to 4 characters), replacing one with the same tag.
    void put(std::string_view tag, const void* data, size_t bytes);
    const std::vector<char>* find(std::string_view tag) const;
    bool empty() const { return m_chunks.empty(); }
    size_t bytes() const;  // in all its chunks

    std::vector<char> serialize() const;
    bool deserialize(const std::vector<char>& bytes, std::string& error);

private:
    std::vector<std::pair<std::string, std::vector<char>>> m_chunks;
};

// The simulation's CPU side as JSON, its history's samples into blob. The bodies'
// pending debris must have been handed to the gas already (the app does so every frame).
json::Value save(const Settings& settings, const Camera& camera, const BodySystem& bodies, Blob& blob);

// Reads that back. Fields a file lacks (one from an older version) keep the values the
// arguments had. On error nothing is changed and error says why.
bool load(const json::Value& root, const Blob& blob, Settings& settings, Camera& camera, BodySystem& bodies,
          std::string& error);

// The settings the simulation reads, as text (the rest only change how it looks: the
// timeline replays changes to these; src/timeline.h). And applied back.
std::string simulationSettings(const Settings& settings);
bool applySimulationSettings(const std::string& text, Settings& settings);

// NAME.bhscene.bin beside NAME.bhscene.
std::filesystem::path sidecarPath(const std::filesystem::path& path);
// Writes both files (each to a temporary first, so a failed save leaves the old ones).
// The JSON records the sidecar's size and checksum.
bool writeFiles(const std::filesystem::path& path, json::Value root, const Blob& blob, std::string& error);
// Reads both and checks that they belong together.
bool readFiles(const std::filesystem::path& path, json::Value& root, Blob& blob, std::string& error);

// 64-bit FNV-1a hash, as 16 hex digits.
std::string checksum(const std::vector<char>& bytes);

}  // namespace scene
