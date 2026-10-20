#pragma once

// The system's Open and Save As dialogs, for scene files (Windows; elsewhere they
// return nothing, and scenes are opened from the command line or by dropping them on
// the window).

#include <filesystem>
#include <optional>
#include <string>

struct GLFWwindow;

namespace dialog {

// The scene file chosen, starting in `folder`; nothing if cancelled.
std::optional<std::filesystem::path> openScene(GLFWwindow* owner, const std::filesystem::path& folder);
// Where to save a scene, offering `name` in `folder`; nothing if cancelled.
std::optional<std::filesystem::path> saveScene(GLFWwindow* owner, const std::filesystem::path& folder,
                                               const std::filesystem::path& name);

}  // namespace dialog
