#pragma once

#include <glad/gl.h>
#include <glm/glm.hpp>

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

// A GLSL program built from files. Supports `#include "file"` relative to the
// shader directory and rebuilds itself when any source file changes on disk.
// A failed rebuild keeps the previous working program. `defines` (whole lines such
// as "#define LENSES 1") go right after the #version line: variants of one source.
class ShaderProgram {
public:
    ShaderProgram() = default;
    ShaderProgram(const ShaderProgram&) = delete;
    ShaderProgram& operator=(const ShaderProgram&) = delete;
    ~ShaderProgram();

    bool build(const std::filesystem::path& dir, const std::string& vertexFile, const std::string& fragmentFile,
               const std::string& defines = {});
    bool buildCompute(const std::filesystem::path& dir, const std::string& computeFile);
    bool reloadIfChanged();  // true if a new program was built

    void use() const { glUseProgram(m_program); }
    bool valid() const { return m_program != 0; }
    GLuint id() const { return m_program; }
    const std::string& error() const { return m_error; }
    const std::string& name() const { return m_fragmentFile; }

    GLint location(const char* name);
    void set(const char* name, int v) { glProgramUniform1i(m_program, location(name), v); }
    void set(const char* name, float v) { glProgramUniform1f(m_program, location(name), v); }
    void set(const char* name, const glm::vec2& v) { glProgramUniform2f(m_program, location(name), v.x, v.y); }
    void set(const char* name, const glm::vec3& v) { glProgramUniform3f(m_program, location(name), v.x, v.y, v.z); }
    void set(const char* name, const glm::ivec2& v) { glProgramUniform2i(m_program, location(name), v.x, v.y); }
    void set(const char* name, const glm::ivec3& v) { glProgramUniform3i(m_program, location(name), v.x, v.y, v.z); }
    void set(const char* name, const glm::vec4& v) {
        glProgramUniform4f(m_program, location(name), v.x, v.y, v.z, v.w);
    }
    void set(const char* name, const glm::mat3& m);
    void set(const char* name, const glm::mat4& m);

private:
    struct Source {
        std::filesystem::path path;
        std::filesystem::file_time_type stamp;
    };
    bool compile();
    bool expand(const std::filesystem::path& file, std::string& out, int depth);

    std::filesystem::path m_dir;
    std::string m_vertexFile, m_fragmentFile;
    std::string m_defines;
    GLuint m_program = 0;
    std::string m_error;
    std::vector<Source> m_sources;  // files that went into the last build attempt
    std::unordered_map<std::string, GLint> m_locations;
};

// A color texture with a framebuffer that renders into it.
struct RenderTarget {
    GLuint fbo = 0;
    GLuint texture = 0;
    int width = 0;
    int height = 0;
    GLenum format = 0;

    void ensure(int w, int h, GLenum internalFormat);  // (re)allocates when the size or format changes
    void destroy();
};
