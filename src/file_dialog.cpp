#include "file_dialog.h"

#include <algorithm>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commdlg.h>
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>
#endif

namespace dialog {

#ifdef _WIN32
namespace {

std::optional<std::filesystem::path> run(GLFWwindow* owner, const std::filesystem::path& folder,
                                         const std::wstring& name, bool save) {
    std::error_code ec;
    std::filesystem::create_directories(folder, ec);
    const std::wstring initial = std::filesystem::absolute(folder, ec).wstring();
    wchar_t file[MAX_PATH * 4] = {};
    name.copy(file, std::min(name.size(), static_cast<size_t>(MAX_PATH * 4 - 1)));

    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = owner ? glfwGetWin32Window(owner) : nullptr;
    ofn.lpstrFilter = L"Black hole scenes (*.bhscene)\0*.bhscene\0All files (*.*)\0*.*\0";
    ofn.lpstrFile = file;
    ofn.nMaxFile = MAX_PATH * 4;
    ofn.lpstrInitialDir = initial.c_str();
    ofn.lpstrDefExt = L"bhscene";
    // The working directory is where the shaders and assets are found: leave it alone.
    ofn.Flags = OFN_NOCHANGEDIR | OFN_PATHMUSTEXIST | (save ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
    const BOOL chosen = save ? GetSaveFileNameW(&ofn) : GetOpenFileNameW(&ofn);
    if (!chosen) return std::nullopt;
    return std::filesystem::path(file);
}

}  // namespace

std::optional<std::filesystem::path> openScene(GLFWwindow* owner, const std::filesystem::path& folder) {
    return run(owner, folder, L"", false);
}

std::optional<std::filesystem::path> saveScene(GLFWwindow* owner, const std::filesystem::path& folder,
                                               const std::filesystem::path& name) {
    return run(owner, folder, name.wstring(), true);
}
#else
std::optional<std::filesystem::path> openScene(GLFWwindow*, const std::filesystem::path&) {
    return std::nullopt;
}

std::optional<std::filesystem::path> saveScene(GLFWwindow*, const std::filesystem::path&, const std::filesystem::path&) {
    return std::nullopt;
}
#endif

}  // namespace dialog
