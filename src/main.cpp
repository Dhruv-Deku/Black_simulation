#include "app.h"
#include "config.h"

#include <cstdio>
#include <filesystem>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

// Ask hybrid-graphics laptops to run us on the discrete GPU.
extern "C" {
__declspec(dllexport) unsigned long NvOptimusEnablement = 1;
__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
}
#endif

namespace fs = std::filesystem;

// Prefers the source tree's shaders (so edits hot-reload while developing), then
// the copy next to the executable.
static fs::path findShaderDir() {
    fs::path source = BH_SHADER_DIR;
    if (fs::exists(source / "trace.frag")) return source;
#ifdef _WIN32
    wchar_t exe[MAX_PATH];
    GetModuleFileNameW(nullptr, exe, MAX_PATH);
    return fs::path(exe).parent_path() / "shaders";
#else
    return fs::current_path() / "shaders";
#endif
}

int main(int argc, char** argv) {
    LaunchOptions options;
    std::string error;
    if (!parseArguments(argc, argv, options, error)) {
        std::fprintf(stderr, "%s\n\n", error.c_str());
        printUsage();
        return 1;
    }
    if (options.showHelp) {
        printUsage();
        return 0;
    }
    App app;
    return app.run(options, findShaderDir());
}
