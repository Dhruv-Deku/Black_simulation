cmake_minimum_required(VERSION 3.24)
project(BlackHoleSimulator LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(FetchContent)

# ---------------------------------------------------------------------------
# Third-party dependencies (downloaded on first configure)
# ---------------------------------------------------------------------------
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
FetchContent_Declare(glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
    GIT_SHALLOW    TRUE)

set(GLM_BUILD_LIBRARY OFF CACHE BOOL "" FORCE)
set(GLM_BUILD_TESTS   OFF CACHE BOOL "" FORCE)
set(GLM_BUILD_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG        1.0.3
    GIT_SHALLOW    TRUE)

FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.92.9b
    GIT_SHALLOW    TRUE)

FetchContent_Declare(stb
    GIT_REPOSITORY https://github.com/nothings/stb.git
    GIT_TAG        2c980bb59875b0d32144a71867fbdebb2f77cd20)

FetchContent_MakeAvailable(glfw glm imgui stb)

# Dear ImGui has no CMake build of its own.
add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)
target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
target_link_libraries(imgui PUBLIC glfw)

add_library(stb INTERFACE)
target_include_directories(stb INTERFACE ${stb_SOURCE_DIR})

# OpenGL 4.6 core loader, generated with glad2 and vendored in external/glad.
add_library(glad STATIC external/glad/src/gl.c)
target_include_directories(glad PUBLIC external/glad/include)

set_target_properties(glfw imgui glad PROPERTIES FOLDER "external")

# ---------------------------------------------------------------------------
# Project
# ---------------------------------------------------------------------------
if(MSVC)
    set(BH_WARNINGS /W4 /permissive- /utf-8 /MP)
else()
    set(BH_WARNINGS -Wall -Wextra)
endif()

# Kerr physics: no graphics dependencies, shared by the app and the tests.
add_library(bh_physics STATIC
    src/kerr.cpp src/kerr.h
    src/central.cpp src/central.h
    src/debris.cpp src/debris.h
    src/trajectory.cpp src/trajectory.h
    src/tidal.cpp src/tidal.h
    src/gravity.cpp src/gravity.h)
target_include_directories(bh_physics PUBLIC src)
target_compile_options(bh_physics PRIVATE ${BH_WARNINGS})

# Shaders are loaded from the source tree while developing so edits hot-reload,
# and from a copy next to the executable otherwise.
set(BH_SHADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/shaders")
configure_file(src/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/generated/config.h @ONLY)

add_executable(blackhole
    src/main.cpp
    src/app.cpp           src/app.h
    src/bodies.cpp        src/bodies.h
    src/blackbody.cpp     src/blackbody.h
    src/camera.cpp        src/camera.h
    src/gl_util.cpp       src/gl_util.h
    src/particles.cpp     src/particles.h
    src/recorder.cpp      src/recorder.h
    src/renderer.cpp      src/renderer.h
    src/ui.cpp            src/ui.h
    src/settings.h)
target_include_directories(blackhole PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_link_libraries(blackhole PRIVATE bh_physics glad glfw imgui glm::glm stb)
target_compile_options(blackhole PRIVATE ${BH_WARNINGS})
if(WIN32)
    # Media Foundation encodes the MP4 videos.
    target_link_libraries(blackhole PRIVATE opengl32 mfplat mfreadwrite mfuuid ole32)
    target_compile_definitions(blackhole PRIVATE NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()
add_custom_command(TARGET blackhole POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/shaders" "$<TARGET_FILE_DIR:blackhole>/shaders"
    COMMENT "Copying shaders next to the executable")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT blackhole)
set_target_properties(blackhole PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
enable_testing()
add_executable(kerr_tests tests/test_kerr.cpp)
target_link_libraries(kerr_tests PRIVATE bh_physics)
target_compile_options(kerr_tests PRIVATE ${BH_WARNINGS})
add_test(NAME kerr_tests COMMAND kerr_tests)
