#pragma once

// A small JSON reader and writer for the scene files (src/scene_file.h).
//
// Numbers keep the text they were read or written with, so a value goes out and comes
// back bit for bit: doubles are written in their shortest exact form (std::to_chars),
// floats in theirs, and each is parsed back with the matching std::from_chars. Objects
// keep their keys in order, so writing the same scene twice gives the same file.
// JSON has no infinities or NaNs; they are written as the strings "inf", "-inf", "nan".

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace json {

class Value {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Value() = default;
    static Value boolean(bool b);
    static Value number(double d);
    static Value number(float f);
    static Value number(long long i);
    static Value string(std::string s);
    static Value array();
    static Value object();
    // A number from its text as written in a file (not checked here).
    static Value numberText(std::string text);

    Type type() const { return m_type; }
    bool isNull() const { return m_type == Type::Null; }
    bool isObject() const { return m_type == Type::Object; }
    bool isArray() const { return m_type == Type::Array; }

    // Scalars. Each returns false (leaving out alone) if the value is of another type or
    // out of range.
    bool get(bool& out) const;
    bool get(double& out) const;
    bool get(float& out) const;
    bool get(long long& out) const;
    bool get(std::string& out) const;
    const std::string& text() const { return m_text; }  // a number's or a string's text

    // Arrays.
    size_t size() const { return m_items.size(); }
    const Value& operator[](size_t i) const { return m_items[i]; }
    Value& operator[](size_t i) { return m_items[i]; }
    void push(Value v) { m_items.push_back(std::move(v)); }

    // Objects: members in the order they were added.
    const Value* find(std::string_view key) const;
    // Adds a member, or replaces the one with that key.
    Value& set(std::string_view key, Value v);
    const std::vector<std::pair<std::string, Value>>& members() const { return m_members; }

    bool operator==(const Value& other) const;

private:
    Type m_type = Type::Null;
    bool m_bool = false;
    std::string m_text;  // a number's digits or a string's characters
    std::vector<Value> m_items;
    std::vector<std::pair<std::string, Value>> m_members;
};

// Two-space indentation; arrays that hold only numbers, booleans or short strings stay on
// one line.
std::string dump(const Value& value);
// Returns false with a message (and the line it was found on) on malformed input.
bool parse(std::string_view text, Value& out, std::string& error);

}  // namespace json
