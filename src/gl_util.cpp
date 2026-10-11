#include "gl_util.h"

#include <glm/gtc/type_ptr.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace {

bool readFile(const fs::path& path, std::string& out) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::ostringstream ss;
    ss << in.rdbuf();
    out = ss.str();
    return true;
}

GLuint compileStage(GLenum stage, const std::string& source, const std::string& label, std::string& error) {
    GLuint shader = glCreateShader(stage);
    const char* text = source.c_str();
    glShaderSource(shader, 1, &text, nullptr);
    glCompileShader(shader);
    GLint ok = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint length = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
        std::string log(static_cast<size_t>(std::max(length, 1)), '\0');
        glGetShaderInfoLog(shader, length, nullptr, log.data());
        error += label + ":\n" + log.c_str() + "\n";
        glDeleteShader(shader);
        return 0;
    }
    return shader;
}

fs::file_time_type stampOf(const fs::path& path) {
    std::error_code ec;
    auto t = fs::last_write_time(path, ec);
    return ec ? fs::file_time_type::min() : t;
}

}  // namespace

ShaderProgram::~ShaderProgram() {
    if (m_program) glDeleteProgram(m_program);
}

bool ShaderProgram::build(const fs::path& dir, const std::string& vertexFile, const std::string& fragmentFile) {
    m_dir = dir;
    m_vertexFile = vertexFile;
    m_fragmentFile = fragmentFile;
    return compile();
}

bool ShaderProgram::buildCompute(const fs::path& dir, const std::string& computeFile) {
    m_dir = dir;
    m_vertexFile.clear();  // an empty vertex stage marks a compute program
    m_fragmentFile = computeFile;
    return compile();
}

bool ShaderProgram::expand(const fs::path& file, std::string& out, int depth) {
    if (depth > 8) {
        m_error += "include depth exceeded at " + file.string() + "\n";
        return false;
    }
    std::string text;
    if (!readFile(file, text)) {
        m_error += "cannot read " + file.string() + "\n";
        return false;
    }
    m_sources.push_back({file, stampOf(file)});

    // Each file gets its own GLSL source-string number so error messages can be
    // traced back: "0(12)" is line 12 of the first file listed in the error.
    const int fileIndex = static_cast<int>(m_sources.size()) - 1;
    std::istringstream lines(text);
    std::string line;
    int lineNo = 0;
    while (std::getline(lines, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        size_t pos = line.find_first_not_of(" \t");
        if (pos != std::string::npos && line.compare(pos, 8, "#include") == 0) {
            size_t open = line.find('"', pos), close = line.rfind('"');
            if (open == std::string::npos || close <= open) {
                m_error += file.filename().string() + ":" + std::to_string(lineNo) + ": malformed #include\n";
                return false;
            }
            if (!expand(m_dir / line.substr(open + 1, close - open - 1), out, depth + 1)) return false;
            out += "#line " + std::to_string(lineNo + 1) + " " + std::to_string(fileIndex) + "\n";
            continue;
        }
        out += line;
        out += '\n';
        if (lineNo == 1 && line.rfind("#version", 0) == 0) {
            out += "#line 2 " + std::to_string(fileIndex) + "\n";
        }
    }
    return true;
}

bool ShaderProgram::compile() {
    m_error.clear();
    m_sources.clear();

    const bool compute = m_vertexFile.empty();
    std::string vsSource, fsSource;
    if (!compute && !expand(m_dir / m_vertexFile, vsSource, 0)) return false;
    const size_t vsFiles = m_sources.size();
    if (!expand(m_dir / m_fragmentFile, fsSource, 0)) return false;

    std::string fileList = "  files:";
    for (size_t i = vsFiles; i < m_sources.size(); ++i) {
        fileList += " [" + std::to_string(i) + "] " + m_sources[i].path.filename().string();
    }

    GLuint vs = compute ? 0 : compileStage(GL_VERTEX_SHADER, vsSource, m_vertexFile, m_error);
    GLuint fs = compileStage(compute ? GL_COMPUTE_SHADER : GL_FRAGMENT_SHADER, fsSource, m_fragmentFile + fileList,
                             m_error);
    if ((!compute && !vs) || !fs) {
        if (vs) glDeleteShader(vs);
        if (fs) glDeleteShader(fs);
        return false;
    }

    GLuint program = glCreateProgram();
    if (vs) glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);
    if (vs) glDeleteShader(vs);
    glDeleteShader(fs);
    GLint ok = GL_FALSE;
    glGetProgramiv(program, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[4096];
        glGetProgramInfoLog(program, sizeof(log), nullptr, log);
        m_error += "link " + m_fragmentFile + ":\n" + log;
        glDeleteProgram(program);
        return false;
    }

    if (m_program) glDeleteProgram(m_program);
    m_program = program;
    m_locations.clear();
    return true;
}

bool ShaderProgram::reloadIfChanged() {
    for (const Source& s : m_sources) {
        if (stampOf(s.path) != s.stamp) return compile();
    }
    return false;
}

GLint ShaderProgram::location(const char* name) {
    auto it = m_locations.find(name);
    if (it != m_locations.end()) return it->second;
    GLint loc = glGetUniformLocation(m_program, name);
    m_locations.emplace(name, loc);
    return loc;
}

void ShaderProgram::set(const char* name, const glm::mat3& m) {
    glProgramUniformMatrix3fv(m_program, location(name), 1, GL_FALSE, glm::value_ptr(m));
}

void RenderTarget::ensure(int w, int h, GLenum internalFormat) {
    w = std::max(w, 1);
    h = std::max(h, 1);
    if (fbo && w == width && h == height && internalFormat == format) return;
    destroy();
    width = w;
    height = h;
    format = internalFormat;

    glCreateTextures(GL_TEXTURE_2D, 1, &texture);
    glTextureStorage2D(texture, 1, internalFormat, w, h);
    glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glCreateFramebuffers(1, &fbo);
    glNamedFramebufferTexture(fbo, GL_COLOR_ATTACHMENT0, texture, 0);
}

void RenderTarget::destroy() {
    if (fbo) glDeleteFramebuffers(1, &fbo);
    if (texture) glDeleteTextures(1, &texture);
    fbo = texture = 0;
    width = height = 0;
}

void ShaderProgram::set(const char* name, const glm::mat4& m) {
    glProgramUniformMatrix4fv(m_program, location(name), 1, GL_FALSE, &m[0][0]);
}
