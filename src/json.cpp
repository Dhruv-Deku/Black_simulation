#include "json.h"

#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <system_error>

namespace json {

namespace {

// Shortest text that reads back as exactly this value.
template <class T>
std::string shortest(T value) {
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    return std::string(buffer, result.ptr);
}

template <class T>
bool readNumber(const std::string& text, T& out) {
    T value{};
    const char* end = text.data() + text.size();
    const auto result = std::from_chars(text.data(), end, value);
    if (result.ec != std::errc() || result.ptr != end) return false;
    out = value;
    return true;
}

// Infinities and NaNs travel as strings.
bool readSpecial(const std::string& text, double& out) {
    if (text == "inf") out = INFINITY;
    else if (text == "-inf") out = -INFINITY;
    else if (text == "nan") out = NAN;
    else return false;
    return true;
}

void appendQuoted(std::string& out, const std::string& s) {
    out += '"';
    for (const char ch : s) {
        const unsigned char c = static_cast<unsigned char>(ch);
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char code[8];
                    std::snprintf(code, sizeof(code), "\\u%04x", c);
                    out += code;
                } else {
                    out += ch;  // UTF-8 passes through
                }
        }
    }
    out += '"';
}

bool isScalar(const Value& v) {
    return v.type() != Value::Type::Array && v.type() != Value::Type::Object &&
           (v.type() != Value::Type::String || v.text().size() < 40);
}

void write(std::string& out, const Value& v, int depth) {
    const std::string pad(static_cast<size_t>(2 * (depth + 1)), ' ');
    const std::string closePad(static_cast<size_t>(2 * depth), ' ');
    switch (v.type()) {
        case Value::Type::Null: out += "null"; return;
        case Value::Type::Bool: {
            bool b = false;
            v.get(b);
            out += b ? "true" : "false";
            return;
        }
        case Value::Type::Number: out += v.text(); return;
        case Value::Type::String: appendQuoted(out, v.text()); return;
        case Value::Type::Array: {
            if (v.size() == 0) {
                out += "[]";
                return;
            }
            bool flat = true;
            for (size_t i = 0; i < v.size() && flat; ++i) flat = isScalar(v[i]);
            if (flat) {
                out += '[';
                for (size_t i = 0; i < v.size(); ++i) {
                    if (i) out += ", ";
                    write(out, v[i], depth + 1);
                }
                out += ']';
                return;
            }
            out += "[\n";
            for (size_t i = 0; i < v.size(); ++i) {
                out += pad;
                write(out, v[i], depth + 1);
                out += i + 1 < v.size() ? ",\n" : "\n";
            }
            out += closePad + "]";
            return;
        }
        case Value::Type::Object: {
            if (v.members().empty()) {
                out += "{}";
                return;
            }
            out += "{\n";
            for (size_t i = 0; i < v.members().size(); ++i) {
                out += pad;
                appendQuoted(out, v.members()[i].first);
                out += ": ";
                write(out, v.members()[i].second, depth + 1);
                out += i + 1 < v.members().size() ? ",\n" : "\n";
            }
            out += closePad + "}";
            return;
        }
    }
}

// Recursive descent over the text, counting lines for the error messages.
class Parser {
public:
    explicit Parser(std::string_view text) : m_text(text) {}

    bool document(Value& out, std::string& error) {
        skipSpace();
        if (!value(out, 0)) {
            error = m_error;
            return false;
        }
        skipSpace();
        if (m_pos != m_text.size()) {
            fail("unexpected text after the end");
            error = m_error;
            return false;
        }
        return true;
    }

private:
    static constexpr int kMaxDepth = 64;

    bool fail(const char* what) {
        if (m_error.empty()) m_error = std::string(what) + " on line " + std::to_string(m_line);
        return false;
    }
    bool atEnd() const { return m_pos >= m_text.size(); }
    char peek() const { return atEnd() ? '\0' : m_text[m_pos]; }
    void skipSpace() {
        while (!atEnd()) {
            const char c = m_text[m_pos];
            if (c == '\n') ++m_line;
            if (c != ' ' && c != '\t' && c != '\n' && c != '\r') break;
            ++m_pos;
        }
    }
    bool literal(std::string_view word) {
        if (m_text.substr(m_pos, word.size()) != word) return fail("unknown word");
        m_pos += word.size();
        return true;
    }

    bool value(Value& out, int depth) {
        if (depth > kMaxDepth) return fail("nested too deeply");
        switch (peek()) {
            case '{': return object(out, depth);
            case '[': return array(out, depth);
            case '"': {
                std::string s;
                if (!string(s)) return false;
                out = Value::string(std::move(s));
                return true;
            }
            case 't': out = Value::boolean(true); return literal("true");
            case 'f': out = Value::boolean(false); return literal("false");
            case 'n': out = Value(); return literal("null");
            default: return number(out);
        }
    }

    bool number(Value& out) {
        const size_t start = m_pos;
        if (peek() == '-') ++m_pos;
        while (!atEnd()) {
            const char c = m_text[m_pos];
            if ((c >= '0' && c <= '9') || c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-') ++m_pos;
            else break;
        }
        std::string text(m_text.substr(start, m_pos - start));
        double check = 0.0;
        if (text.empty() || !readNumber(text, check)) {
            // from_chars rejects out-of-range values; JSON allows them, but nothing here writes them.
            return fail(text.empty() ? "expected a value" : "malformed number");
        }
        out = Value::numberText(std::move(text));
        return true;
    }

    static void appendUtf8(std::string& out, uint32_t code) {
        if (code < 0x80) {
            out += static_cast<char>(code);
        } else if (code < 0x800) {
            out += static_cast<char>(0xC0 | (code >> 6));
            out += static_cast<char>(0x80 | (code & 0x3F));
        } else if (code < 0x10000) {
            out += static_cast<char>(0xE0 | (code >> 12));
            out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (code & 0x3F));
        } else {
            out += static_cast<char>(0xF0 | (code >> 18));
            out += static_cast<char>(0x80 | ((code >> 12) & 0x3F));
            out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (code & 0x3F));
        }
    }

    bool hex4(uint32_t& code) {
        if (m_pos + 4 > m_text.size()) return fail("short \\u escape");
        code = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = m_text[m_pos++];
            code <<= 4;
            if (c >= '0' && c <= '9') code |= static_cast<uint32_t>(c - '0');
            else if (c >= 'a' && c <= 'f') code |= static_cast<uint32_t>(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') code |= static_cast<uint32_t>(c - 'A' + 10);
            else return fail("bad \\u escape");
        }
        return true;
    }

    bool string(std::string& out) {
        ++m_pos;  // opening quote
        while (true) {
            if (atEnd()) return fail("unterminated string");
            const char c = m_text[m_pos++];
            if (c == '"') return true;
            if (c == '\n') return fail("line break inside a string");
            if (c != '\\') {
                out += c;
                continue;
            }
            if (atEnd()) return fail("unterminated string");
            const char e = m_text[m_pos++];
            switch (e) {
                case '"': out += '"'; break;
                case '\\': out += '\\'; break;
                case '/': out += '/'; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'n': out += '\n'; break;
                case 'r': out += '\r'; break;
                case 't': out += '\t'; break;
                case 'u': {
                    uint32_t code = 0;
                    if (!hex4(code)) return false;
                    // A UTF-16 surrogate pair spells one character beyond the first 65,536.
                    if (code >= 0xD800 && code < 0xDC00 && m_text.substr(m_pos, 2) == "\\u") {
                        m_pos += 2;
                        uint32_t low = 0;
                        if (!hex4(low)) return false;
                        code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
                    }
                    appendUtf8(out, code);
                    break;
                }
                default: return fail("bad escape in a string");
            }
        }
    }

    bool array(Value& out, int depth) {
        ++m_pos;
        out = Value::array();
        skipSpace();
        if (peek() == ']') {
            ++m_pos;
            return true;
        }
        while (true) {
            skipSpace();
            Value item;
            if (!value(item, depth + 1)) return false;
            out.push(std::move(item));
            skipSpace();
            if (peek() == ',') {
                ++m_pos;
                continue;
            }
            if (peek() == ']') {
                ++m_pos;
                return true;
            }
            return fail("expected , or ] in an array");
        }
    }

    bool object(Value& out, int depth) {
        ++m_pos;
        out = Value::object();
        skipSpace();
        if (peek() == '}') {
            ++m_pos;
            return true;
        }
        while (true) {
            skipSpace();
            if (peek() != '"') return fail("expected a key in quotes");
            std::string key;
            if (!string(key)) return false;
            skipSpace();
            if (peek() != ':') return fail("expected : after a key");
            ++m_pos;
            skipSpace();
            Value item;
            if (!value(item, depth + 1)) return false;
            out.set(key, std::move(item));
            skipSpace();
            if (peek() == ',') {
                ++m_pos;
                continue;
            }
            if (peek() == '}') {
                ++m_pos;
                return true;
            }
            return fail("expected , or } in an object");
        }
    }

    std::string_view m_text;
    size_t m_pos = 0;
    int m_line = 1;
    std::string m_error;
};

}  // namespace

Value Value::boolean(bool b) {
    Value v;
    v.m_type = Type::Bool;
    v.m_bool = b;
    return v;
}

Value Value::number(double d) {
    if (!std::isfinite(d)) return string(std::isnan(d) ? "nan" : (d > 0.0 ? "inf" : "-inf"));
    return numberText(shortest(d));
}

Value Value::number(float f) {
    if (!std::isfinite(f)) return number(static_cast<double>(f));
    return numberText(shortest(f));
}

Value Value::number(long long i) {
    return numberText(shortest(i));
}

Value Value::string(std::string s) {
    Value v;
    v.m_type = Type::String;
    v.m_text = std::move(s);
    return v;
}

Value Value::array() {
    Value v;
    v.m_type = Type::Array;
    return v;
}

Value Value::object() {
    Value v;
    v.m_type = Type::Object;
    return v;
}

Value Value::numberText(std::string text) {
    Value v;
    v.m_type = Type::Number;
    v.m_text = std::move(text);
    return v;
}

bool Value::get(bool& out) const {
    if (m_type != Type::Bool) return false;
    out = m_bool;
    return true;
}

bool Value::get(double& out) const {
    if (m_type == Type::String) return readSpecial(m_text, out);
    return m_type == Type::Number && readNumber(m_text, out);
}

bool Value::get(float& out) const {
    if (m_type == Type::String) {
        double d = 0.0;
        if (!readSpecial(m_text, d)) return false;
        out = static_cast<float>(d);
        return true;
    }
    return m_type == Type::Number && readNumber(m_text, out);
}

bool Value::get(long long& out) const {
    if (m_type != Type::Number) return false;
    if (readNumber(m_text, out)) return true;
    // Written as a floating-point number (by hand, say): accept it if it is whole.
    double d = 0.0;
    if (!readNumber(m_text, d) || d != std::floor(d) || std::abs(d) > 9.0e18) return false;
    out = static_cast<long long>(d);
    return true;
}

bool Value::get(std::string& out) const {
    if (m_type != Type::String) return false;
    out = m_text;
    return true;
}

const Value* Value::find(std::string_view key) const {
    for (const auto& [k, v] : m_members) {
        if (k == key) return &v;
    }
    return nullptr;
}

Value& Value::set(std::string_view key, Value v) {
    m_type = Type::Object;
    for (auto& [k, existing] : m_members) {
        if (k == key) {
            existing = std::move(v);
            return existing;
        }
    }
    m_members.emplace_back(std::string(key), std::move(v));
    return m_members.back().second;
}

bool Value::operator==(const Value& other) const {
    return m_type == other.m_type && m_bool == other.m_bool && m_text == other.m_text && m_items == other.m_items &&
           m_members == other.m_members;
}

std::string dump(const Value& value) {
    std::string out;
    write(out, value, 0);
    out += '\n';
    return out;
}

bool parse(std::string_view text, Value& out, std::string& error) {
    // A UTF-8 byte order mark (Notepad adds one) is not JSON; skip it.
    if (text.substr(0, 3) == "\xEF\xBB\xBF") text.remove_prefix(3);
    Parser parser(text);
    return parser.document(out, error);
}

}  // namespace json
