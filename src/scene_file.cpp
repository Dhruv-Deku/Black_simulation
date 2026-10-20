#include "scene_file.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <type_traits>

namespace scene {

// What a saved scene must carry that its classes keep private: the bodies' random
// number generator (it places the next preset), their name counter and smoothed frame
// time (it sets the history's step), and each fallback's rate table.
struct Access {
    static std::mt19937& rng(BodySystem& b) { return b.m_rng; }
    static int& counter(BodySystem& b) { return b.m_counter; }
    static double& frameTime(BodySystem& b) { return b.m_frameTime; }
    static double& logT0(tde::Fallback& f) { return f.m_logT0; }
    static double& dLogT(tde::Fallback& f) { return f.m_dLogT; }
    static std::vector<double>& rate(tde::Fallback& f) { return f.m_rate; }
    static std::vector<double>& cumulative(tde::Fallback& f) { return f.m_cumulative; }
};

namespace {

using json::Value;

constexpr char kFormat[] = "bhscene";
constexpr char kBlobMagic[8] = {'B', 'H', 'S', 'C', 'N', 'B', 'I', 'N'};
constexpr uint32_t kBlobVersion = 1;
// One slot of the light-delay history: its kRows rows of kColumns samples.
constexpr size_t kHistorySlotBytes = sizeof(glm::vec4) * History::kRows * History::kColumns;

// Enumerations are stored as numbers; a number past the last enumerator is refused.
template <class E>
constexpr long long kEnumCount = 0;  // unchecked
template <>
constexpr long long kEnumCount<BodyKind> = static_cast<long long>(BodyKind::Moon) + 1;
template <>
constexpr long long kEnumCount<Palette> = static_cast<long long>(Palette::Scorched) + 1;
template <>
constexpr long long kEnumCount<central::Kind> = static_cast<long long>(central::Kind::GasGiant) + 1;
template <>
constexpr long long kEnumCount<central::Host> = static_cast<long long>(central::Host::GasGiant) + 1;
template <>
constexpr long long kEnumCount<CameraMotion> = static_cast<long long>(CameraMotion::RideAlong) + 1;

// Where in the file a bad value was found ("bodies.state.r"), built up as decoding unwinds.
std::string g_where;

// --- Single values -----------------------------------------------------------------
// Every type has an encode and a decode, all declared here first so that the templates
// below find each other.

Value encode(bool v);
Value encode(int v);
Value encode(long long v);
Value encode(float v);
Value encode(double v);
Value encode(const std::string& v);
template <class E>
    requires std::is_enum_v<E>
Value encode(E v);
Value encode(const glm::vec3& v);
template <class T, size_t N>
Value encode(const std::array<T, N>& v);
template <class T>
Value encode(const std::vector<T>& v);
template <class T>
Value encode(const std::deque<T>& v);
Value encode(const central::Object& v);
Value encode(const kerr::Constants& v);
Value encode(const kerr::State& v);
Value encode(const tidal::Star& v);
Value encode(const Body& v);
Value encode(const Impact& v);
Value encode(const BodySystem::Departed& v);
Value encode(const tde::Fallback& v);
Value encode(const tde::Event& v);
Value encode(const penrose::Split& v);

bool decode(const Value& x, bool& v);
bool decode(const Value& x, int& v);
bool decode(const Value& x, long long& v);
bool decode(const Value& x, float& v);
bool decode(const Value& x, double& v);
bool decode(const Value& x, std::string& v);
template <class E>
    requires std::is_enum_v<E>
bool decode(const Value& x, E& v);
bool decode(const Value& x, glm::vec3& v);
template <class T, size_t N>
bool decode(const Value& x, std::array<T, N>& v);
template <class T>
bool decode(const Value& x, std::vector<T>& v);
template <class T>
bool decode(const Value& x, std::deque<T>& v);
bool decode(const Value& x, central::Object& v);
bool decode(const Value& x, kerr::Constants& v);
bool decode(const Value& x, kerr::State& v);
bool decode(const Value& x, tidal::Star& v);
bool decode(const Value& x, Body& v);
bool decode(const Value& x, Impact& v);
bool decode(const Value& x, BodySystem::Departed& v);
bool decode(const Value& x, tde::Fallback& v);
bool decode(const Value& x, tde::Event& v);
bool decode(const Value& x, penrose::Split& v);

// One list of fields serves both directions (the fields() functions below): Writer puts
// each into a JSON object, Reader takes each back if the object has it.
class Writer {
public:
    static constexpr bool kReading = false;
    Value object = Value::object();
    template <class T>
    void operator()(const char* key, const T& v) {
        object.set(key, encode(v));
    }
};

class Reader {
public:
    static constexpr bool kReading = true;
    explicit Reader(const Value& object) : m_object(object) {}
    template <class T>
    void operator()(const char* key, T& v) {
        if (!m_ok) return;
        const Value* x = m_object.find(key);
        if (!x) return;  // not in the file (an older version): keep the default
        if (!decode(*x, v)) {
            g_where = std::string(key) + (g_where.empty() || g_where[0] == '[' ? "" : ".") + g_where;
            m_ok = false;
        }
    }
    bool ok() const { return m_ok; }

private:
    const Value& m_object;
    bool m_ok = true;
};

template <class T>
Value encodeFields(const T& v) {
    Writer w;
    fields(w, const_cast<T&>(v));
    return std::move(w.object);
}

template <class T>
bool decodeFields(const Value& x, T& v) {
    if (!x.isObject()) return false;
    Reader r(x);
    fields(r, v);
    return r.ok();
}

Value encode(bool v) { return Value::boolean(v); }
Value encode(int v) { return Value::number(static_cast<long long>(v)); }
Value encode(long long v) { return Value::number(v); }
Value encode(float v) { return Value::number(v); }
Value encode(double v) { return Value::number(v); }
Value encode(const std::string& v) { return Value::string(v); }
template <class E>
    requires std::is_enum_v<E>
Value encode(E v) {
    return Value::number(static_cast<long long>(v));
}
Value encode(const glm::vec3& v) {
    Value a = Value::array();
    for (int i = 0; i < 3; ++i) a.push(encode(v[i]));
    return a;
}
template <class T, size_t N>
Value encode(const std::array<T, N>& v) {
    Value a = Value::array();
    for (const T& item : v) a.push(encode(item));
    return a;
}
template <class T>
Value encode(const std::vector<T>& v) {
    Value a = Value::array();
    for (const T& item : v) a.push(encode(item));
    return a;
}
template <class T>
Value encode(const std::deque<T>& v) {
    Value a = Value::array();
    for (const T& item : v) a.push(encode(item));
    return a;
}

bool decode(const Value& x, bool& v) { return x.get(v); }
bool decode(const Value& x, int& v) {
    long long i = 0;
    if (!x.get(i) || i < INT_MIN || i > INT_MAX) return false;
    v = static_cast<int>(i);
    return true;
}
bool decode(const Value& x, long long& v) { return x.get(v); }
bool decode(const Value& x, float& v) { return x.get(v); }
bool decode(const Value& x, double& v) { return x.get(v); }
bool decode(const Value& x, std::string& v) { return x.get(v); }
template <class E>
    requires std::is_enum_v<E>
bool decode(const Value& x, E& v) {
    long long i = 0;
    if (!x.get(i) || i < 0 || (kEnumCount<E> > 0 && i >= kEnumCount<E>)) return false;
    v = static_cast<E>(i);
    return true;
}
bool decode(const Value& x, glm::vec3& v) {
    if (!x.isArray() || x.size() != 3) return false;
    for (int i = 0; i < 3; ++i) {
        if (!decode(x[static_cast<size_t>(i)], v[i])) return false;
    }
    return true;
}
template <class T, size_t N>
bool decode(const Value& x, std::array<T, N>& v) {
    if (!x.isArray() || x.size() != N) return false;
    for (size_t i = 0; i < N; ++i) {
        if (!decode(x[i], v[i])) return false;
    }
    return true;
}
template <class T>
bool decode(const Value& x, std::vector<T>& v) {
    if (!x.isArray()) return false;
    std::vector<T> out(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        if (!decode(x[i], out[i])) {
            g_where = "[" + std::to_string(i) + "]" + (g_where.empty() || g_where[0] == '[' ? "" : ".") + g_where;
            return false;
        }
    }
    v = std::move(out);
    return true;
}
template <class T>
bool decode(const Value& x, std::deque<T>& v) {
    std::vector<T> items;
    if (!decode(x, items)) return false;
    v.assign(items.begin(), items.end());
    return true;
}

// --- Fields --------------------------------------------------------------------------

template <class A>
void fields(A& ar, central::Object& o) {
    ar("kind", o.kind);
    ar("name", o.name);
    ar("massSolar", o.massSolar);
    ar("radius", o.radius);
    ar("chi", o.chi);
    ar("inertia", o.inertia);
    ar("surfaceTemperature", o.surfaceTemperature);
    ar("magneticField", o.magneticField);
    ar("magneticTiltDeg", o.magneticTiltDeg);
    ar("pulsar", o.pulsar);
    ar("magnetar", o.magnetar);
    ar("beacon", o.beacon);
    ar("naked", o.naked);
}

template <class A>
void fields(A& ar, kerr::Constants& c) {
    ar("E", c.E);
    ar("L", c.L);
    ar("Q", c.Q);
    ar("mu", c.mu);
}

template <class A>
void fields(A& ar, kerr::State& s) {
    ar("t", s.t);
    ar("r", s.r);
    ar("theta", s.theta);
    ar("phi", s.phi);
    ar("vr", s.vr);
    ar("vtheta", s.vtheta);
}

template <class A>
void fields(A& ar, tidal::Star& s) {
    ar("F", s.F);
    ar("Fdot", s.Fdot);
}

template <class A>
void fields(A& ar, Body& b) {
    ar("name", b.name);
    ar("kind", b.kind);
    ar("constants", b.constants);
    ar("state", b.state);
    ar("mass", b.mass);
    ar("radius", b.radius);
    ar("temperature", b.temperature);
    ar("interiorTemperature", b.interiorTemperature);
    ar("brightness", b.brightness);
    ar("seed", b.seed);
    ar("palette", b.palette);
    ar("spinRate", b.spinRate);
    ar("tidalRadius", b.tidalRadius);
    ar("shape", b.shape);
    ar("maxStretch", b.maxStretch);
    ar("inside", b.inside);
    ar("slot", b.slot);
    ar("trail", b.trail);
    ar("splitRadius", b.splitRadius);
    ar("penroseSplit", b.penroseSplit);
}

template <class A>
void fields(A& ar, penrose::Split& s) {
    ar("valid", s.valid);
    ar("why", s.why);
    ar("r", s.r);
    ar("theta", s.theta);
    ar("m0", s.m0);
    ar("m1", s.m1);
    ar("m2", s.m2);
    ar("E0", s.E0);
    ar("E1", s.E1);
    ar("E2", s.E2);
    ar("L0", s.L0);
    ar("L1", s.L1);
    ar("L2", s.L2);
    std::array<double, 3> v1, v2;
    std::copy(std::begin(s.v1), std::end(s.v1), v1.begin());
    std::copy(std::begin(s.v2), std::end(s.v2), v2.begin());
    ar("v1", v1);
    ar("v2", v2);
    if constexpr (A::kReading) {
        std::copy(v1.begin(), v1.end(), std::begin(s.v1));
        std::copy(v2.begin(), v2.end(), std::begin(s.v2));
    }
    ar("efficiency", s.efficiency);
    ar("spinBefore", s.spinBefore);
    ar("absorbed", s.absorbed);
    ar("applied", s.applied);
}

template <class A>
void fields(A& ar, Impact& i) {
    ar("theta", i.theta);
    ar("phi", i.phi);
    ar("time", i.time);
    ar("energy", i.energy);
    ar("mass", i.mass);
}

template <class A>
void fields(A& ar, BodySystem::Departed& d) {
    ar("body", d.body);
    ar("lostTime", d.lostTime);
}

template <class A>
void fields(A& ar, tde::Fallback& f) {
    ar("boundFraction", f.boundFraction);
    ar("firstReturn", f.firstReturn);
    ar("peakTime", f.peakTime);
    ar("peakRate", f.peakRate);
    ar("lateSlope", f.lateSlope);
    ar("logT0", Access::logT0(f));
    ar("dLogT", Access::dLogT(f));
    ar("rate", Access::rate(f));
    ar("cumulative", Access::cumulative(f));
}

template <class A>
void fields(A& ar, tde::Event& e) {
    ar("time", e.time);
    ar("axis", e.axis);
    ar("pericenter", e.pericenter);
    ar("circularization", e.circularization);
    ar("starMass", e.starMass);
    ar("starRadius", e.starRadius);
    ar("tidalRadius", e.tidalRadius);
    ar("unboundFraction", e.unboundFraction);
    ar("plungingFraction", e.plungingFraction);
    ar("fallback", e.fallback);
}

template <class A>
void fields(A& ar, Settings& s) {
    ar("central", s.central);
    ar("physicalUnits", s.physicalUnits);
    ar("radiusFollowsMass", s.radiusFollowsMass);
    ar("autoStage", s.autoStage);
    ar("accretionRate", s.accretionRate);
    ar("magnetosphere", s.magnetosphere);
    ar("fieldLines", s.fieldLines);
    ar("pulsarBeams", s.pulsarBeams);
    ar("starquakes", s.starquakes);
    ar("magnetosphereGlow", s.magnetosphereGlow);
    ar("surfaceDetail", s.surfaceDetail);
    ar("surfaceBrightness", s.surfaceBrightness);
    ar("corona", s.corona);
    ar("coronaBoost", s.coronaBoost);
    ar("companionStar", s.companionStar);
    ar("filmColor", s.filmColor);
    ar("planetMaps", s.planetMaps);
    ar("spin", s.spin);
    ar("diskEnabled", s.diskEnabled);
    ar("diskPrograde", s.diskPrograde);
    ar("diskInnerAtIsco", s.diskInnerAtIsco);
    ar("diskInner", s.diskInner);
    ar("diskOuter", s.diskOuter);
    ar("diskTemperature", s.diskTemperature);
    ar("diskPhysicalTemperature", s.diskPhysicalTemperature);
    ar("diskBrightness", s.diskBrightness);
    ar("diskOpacity", s.diskOpacity);
    ar("diskTurbulence", s.diskTurbulence);
    ar("relativity", s.relativity);
    ar("lightDelay", s.lightDelay);
    ar("bodyGravity", s.bodyGravity);
    ar("debrisBrightness", s.debrisBrightness);
    ar("debrisOpacity", s.debrisOpacity);
    ar("gasPhysics", s.gasPhysics);
    ar("gasCollisions", s.gasCollisions);
    ar("gasCoolingTime", s.gasCoolingTime);
    ar("gasFloorTemperature", s.gasFloorTemperature);
    ar("debrisGlowFloor", s.debrisGlowFloor);
    ar("tdeFlare", s.tdeFlare);
    ar("outflowFraction", s.outflowFraction);
    ar("outflowSpeed", s.outflowSpeed);
    ar("outflowOpeningDeg", s.outflowOpeningDeg);
    ar("outflowWidening", s.outflowWidening);
    ar("outflowThickness", s.outflowThickness);
    ar("outflowBrightness", s.outflowBrightness);
    ar("flareDisk", s.flareDisk);
    ar("flareCinematic", s.flareCinematic);
    ar("debrisFilaments", s.debrisFilaments);
    ar("flareFastForward", s.flareFastForward);
    ar("flareCamera", s.flareCamera);
    ar("jetsEnabled", s.jetsEnabled);
    ar("jetPower", s.jetPower);
    ar("jetSpeed", s.jetSpeed);
    ar("jetOpeningDeg", s.jetOpeningDeg);
    ar("showPaths", s.showPaths);
    ar("showPredictions", s.showPredictions);
    ar("predictionTime", s.predictionTime);
    ar("pathWidth", s.pathWidth);
    ar("skyMode", s.skyMode);
    ar("starBrightness", s.starBrightness);
    ar("galaxyBrightness", s.galaxyBrightness);
    ar("skyImageBrightness", s.skyImageBrightness);
    ar("renderScale", s.renderScale);
    ar("maxSteps", s.maxSteps);
    ar("stepSize", s.stepSize);
    ar("escapeRadius", s.escapeRadius);
    ar("exposure", s.exposure);
    ar("bloomStrength", s.bloomStrength);
    ar("tonemapper", s.tonemapper);
    ar("debugView", s.debugView);
    ar("upscaler", s.upscaler);
    ar("sharpness", s.sharpness);
    ar("dynamicResolution", s.dynamicResolution);
    ar("targetMs", s.targetMs);
    ar("timeScale", s.timeScale);
    ar("penroseSpinDown", s.penroseSpinDown);
    ar("paused", s.paused);
    ar("sound", s.sound);
    ar("soundVolume", s.soundVolume);
    ar("soundBodies", s.soundBodies);
    ar("soundPulsars", s.soundPulsars);
    ar("soundRumble", s.soundRumble);
    ar("soundWaves", s.soundWaves);
    ar("rayInspector", s.rayInspector);
    ar("slowNearHole", s.slowNearHole);
    ar("autoOrbit", s.autoOrbit);
    ar("autoOrbitSpeed", s.autoOrbitSpeed);
    // Not vsync: that belongs to the screen, not the scene.
}

// What the simulation reads from the settings: the central object, what sets off flashes
// and flares, the bodies' gravity and the gas physics. The time scale is not here: the
// timeline logs each step's length itself.
template <class A>
void simulationFields(A& ar, Settings& s) {
    ar("central", s.central);
    ar("physicalUnits", s.physicalUnits);
    ar("radiusFollowsMass", s.radiusFollowsMass);
    ar("autoStage", s.autoStage);
    ar("accretionRate", s.accretionRate);
    ar("starquakes", s.starquakes);
    ar("surfaceDetail", s.surfaceDetail);
    ar("diskEnabled", s.diskEnabled);
    ar("lightDelay", s.lightDelay);
    ar("bodyGravity", s.bodyGravity);
    ar("gasPhysics", s.gasPhysics);
    ar("gasCollisions", s.gasCollisions);
    ar("gasCoolingTime", s.gasCoolingTime);
    ar("gasFloorTemperature", s.gasFloorTemperature);
    ar("penroseSpinDown", s.penroseSpinDown);
}

template <class A>
void fields(A& ar, Camera& c) {
    ar("r", c.r);
    ar("theta", c.theta);
    ar("phi", c.phi);
    ar("yaw", c.yaw);
    ar("pitch", c.pitch);
    ar("roll", c.roll);
    ar("fovDeg", c.fovDeg);
    ar("motion", c.motion);
    ar("geodesic", c.geodesic);
    ar("path", c.path);
    ar("properTime", c.properTime);
    ar("ended", c.ended);
    ar("rideSlot", c.rideSlot);
    ar("lookAhead", c.lookAhead);
    ar("surface", c.surface);
    ar("scale", c.scale);
    ar("landed", c.landed);
    ar("landingSpeed", c.landingSpeed);
    ar("drift", c.drift);
}

Value encode(const central::Object& v) { return encodeFields(v); }
Value encode(const kerr::Constants& v) { return encodeFields(v); }
Value encode(const kerr::State& v) { return encodeFields(v); }
Value encode(const tidal::Star& v) { return encodeFields(v); }
Value encode(const Body& v) { return encodeFields(v); }
Value encode(const Impact& v) { return encodeFields(v); }
Value encode(const BodySystem::Departed& v) { return encodeFields(v); }
Value encode(const tde::Fallback& v) { return encodeFields(v); }
Value encode(const tde::Event& v) { return encodeFields(v); }
Value encode(const penrose::Split& v) { return encodeFields(v); }

bool decode(const Value& x, central::Object& v) { return decodeFields(x, v); }
bool decode(const Value& x, kerr::Constants& v) { return decodeFields(x, v); }
bool decode(const Value& x, kerr::State& v) { return decodeFields(x, v); }
bool decode(const Value& x, tidal::Star& v) { return decodeFields(x, v); }
bool decode(const Value& x, Body& v) { return decodeFields(x, v); }
bool decode(const Value& x, Impact& v) { return decodeFields(x, v); }
bool decode(const Value& x, BodySystem::Departed& v) { return decodeFields(x, v); }
bool decode(const Value& x, tde::Fallback& v) { return decodeFields(x, v); }
bool decode(const Value& x, tde::Event& v) { return decodeFields(x, v); }
bool decode(const Value& x, penrose::Split& v) { return decodeFields(x, v); }

// The bodies, everything about them but the history's samples (those go in the blob).
template <class A>
void fields(A& ar, BodySystem& b) {
    ar("time", b.time);
    ar("surface", b.surface);
    ar("surfaceOmega", b.surfaceOmega);
    ar("scale", b.scale);
    ar("massSolar", b.massSolar);
    ar("host", b.host);
    ar("centralMassSolar", b.centralMassSolar);
    ar("centralRadiusKm", b.centralRadiusKm);
    ar("centralTemperature", b.centralTemperature);
    ar("centralSpinPeriod", b.centralSpinPeriod);
    ar("accreted", b.accreted);
    ar("impactCount", b.impactCount);
    ar("impacts", b.impacts);
    ar("bodies", b.bodies);
    ar("departed", b.departed);
    ar("events", b.events);
    ar("tidalEvents", b.tidalEvents);
    ar("penroseSplits", b.penroseSplits);
    ar("historyStep", b.history.step);
    ar("historyNewest", b.history.newest);
    std::array<int, History::kSlots> lengths;
    std::copy(std::begin(b.history.length), std::end(b.history.length), lengths.begin());
    ar("historyLength", lengths);
    ar("nameCounter", Access::counter(b));
    ar("frameTime", Access::frameTime(b));
    // The generator's whole state (624 words and a position), as the standard library writes it.
    std::ostringstream rngOut;
    rngOut << Access::rng(b);
    std::string rng = rngOut.str();
    ar("random", rng);
    if constexpr (A::kReading) {
        std::copy(lengths.begin(), lengths.end(), std::begin(b.history.length));
        std::istringstream rngIn(rng);
        rngIn >> Access::rng(b);
    }
}

// Every body's slot must be a distinct row of the history.
bool checkSlots(const BodySystem& b, std::string& error) {
    bool used[History::kSlots] = {};
    auto claim = [&](const Body& body) {
        if (body.slot < 0 || body.slot >= History::kSlots || used[body.slot]) return false;
        used[body.slot] = true;
        return true;
    };
    for (const Body& body : b.bodies) {
        if (!claim(body)) {
            error = "body \"" + body.name + "\" has a bad or repeated history slot";
            return false;
        }
    }
    for (const BodySystem::Departed& d : b.departed) {
        if (!claim(d.body)) {
            error = "departed body \"" + d.body.name + "\" has a bad or repeated history slot";
            return false;
        }
    }
    for (int slot = 0; slot < History::kSlots; ++slot) {
        if (b.history.length[slot] < 0 || b.history.length[slot] > History::kColumns) {
            error = "bad history length";
            return false;
        }
    }
    if (!(b.history.step > 0.0) || b.bodies.size() > static_cast<size_t>(History::kSlots)) {
        error = "bad history step or too many bodies";
        return false;
    }
    return true;
}

template <class T>
void append(std::vector<char>& out, const T& value) {
    const char* p = reinterpret_cast<const char*>(&value);
    out.insert(out.end(), p, p + sizeof(T));
}

template <class T>
bool take(const std::vector<char>& in, size_t& offset, T& value) {
    if (offset > in.size() || in.size() - offset < sizeof(T)) return false;
    std::memcpy(&value, in.data() + offset, sizeof(T));
    offset += sizeof(T);
    return true;
}

std::string utf8(const std::filesystem::path& path) {
    const std::u8string s = path.u8string();
    return std::string(s.begin(), s.end());
}

bool readWholeFile(const std::filesystem::path& path, std::vector<char>& out) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    in.seekg(0, std::ios::end);
    const std::streamoff size = in.tellg();
    if (size < 0) return false;
    in.seekg(0, std::ios::beg);
    out.resize(static_cast<size_t>(size));
    return size == 0 || static_cast<bool>(in.read(out.data(), size));
}

bool writeWholeFile(const std::filesystem::path& path, const char* data, size_t bytes) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out.write(data, static_cast<std::streamsize>(bytes));
    out.close();
    return static_cast<bool>(out);
}

}  // namespace

// --- Blob ------------------------------------------------------------------------------

void Blob::put(std::string_view tag, const void* data, size_t bytes) {
    const char* p = static_cast<const char*>(data);
    std::vector<char> chunk(p, p + bytes);
    const std::string key(tag.substr(0, 4));
    for (auto& [k, v] : m_chunks) {
        if (k == key) {
            v = std::move(chunk);
            return;
        }
    }
    m_chunks.emplace_back(key, std::move(chunk));
}

const std::vector<char>* Blob::find(std::string_view tag) const {
    for (const auto& [k, v] : m_chunks) {
        if (k == tag) return &v;
    }
    return nullptr;
}

// "BHSCNBIN", a version and a chunk count (uint32 each), then per chunk a 4-byte tag
// (padded with zeros), its length (uint64) and its bytes. Little-endian, as on every
// machine this runs on.
std::vector<char> Blob::serialize() const {
    std::vector<char> out(kBlobMagic, kBlobMagic + sizeof(kBlobMagic));
    append(out, kBlobVersion);
    append(out, static_cast<uint32_t>(m_chunks.size()));
    for (const auto& [tag, bytes] : m_chunks) {
        char name[4] = {};
        std::memcpy(name, tag.data(), std::min<size_t>(tag.size(), 4));
        out.insert(out.end(), name, name + 4);
        append(out, static_cast<uint64_t>(bytes.size()));
        out.insert(out.end(), bytes.begin(), bytes.end());
    }
    return out;
}

bool Blob::deserialize(const std::vector<char>& bytes, std::string& error) {
    m_chunks.clear();
    size_t offset = sizeof(kBlobMagic);
    uint32_t version = 0, count = 0;
    if (bytes.size() < offset || std::memcmp(bytes.data(), kBlobMagic, sizeof(kBlobMagic)) != 0 ||
        !take(bytes, offset, version) || !take(bytes, offset, count)) {
        error = "the sidecar is not a scene's binary file";
        return false;
    }
    if (version > kBlobVersion) {
        error = "the sidecar was written by a newer version of the simulator";
        return false;
    }
    for (uint32_t i = 0; i < count; ++i) {
        if (bytes.size() - offset < 4) {
            error = "the sidecar is cut short";
            return false;
        }
        std::string tag(bytes.data() + offset, 4);
        tag.resize(std::strlen(tag.c_str()));
        offset += 4;
        uint64_t size = 0;
        if (!take(bytes, offset, size) || size > bytes.size() - offset) {
            error = "the sidecar is cut short";
            return false;
        }
        m_chunks.emplace_back(tag, std::vector<char>(bytes.begin() + static_cast<std::ptrdiff_t>(offset),
                                                     bytes.begin() + static_cast<std::ptrdiff_t>(offset + size)));
        offset += static_cast<size_t>(size);
    }
    return true;
}

// --- Scenes ----------------------------------------------------------------------------

json::Value save(const Settings& settings, const Camera& camera, const BodySystem& bodies, Blob& blob) {
    Value root = Value::object();
    root.set("format", Value::string(kFormat));
    root.set("version", Value::number(static_cast<long long>(kVersion)));
    root.set("settings", encodeFields(settings));
    root.set("camera", encodeFields(camera));
    root.set("bodies", encodeFields(bodies));

    // The history rows of every slot in use, in slot order: a slot number (int32) and its
    // rows. Unused slots are left out; the renderer reads only the samples counted in
    // historyLength anyway.
    std::vector<char> history;
    for (int slot = 0; slot < History::kSlots; ++slot) {
        if (bodies.history.length[slot] <= 0) continue;
        append(history, static_cast<int32_t>(slot));
        const char* rows = reinterpret_cast<const char*>(
            &bodies.history.image[static_cast<size_t>(History::kRows * slot) * History::kColumns]);
        history.insert(history.end(), rows, rows + kHistorySlotBytes);
    }
    blob.put("HIST", history.data(), history.size());
    return root;
}

bool load(const json::Value& root, const Blob& blob, Settings& settings, Camera& camera, BodySystem& bodies,
          std::string& error) {
    std::string format;
    const Value* f = root.find("format");
    if (!root.isObject() || !f || !f->get(format) || format != kFormat) {
        error = "not a scene file";
        return false;
    }
    long long version = 0;
    const Value* v = root.find("version");
    if (!v || !v->get(version) || version < 1) {
        error = "the scene file has no version";
        return false;
    }
    if (version > kVersion) {
        error = "the scene was saved by a newer version of the simulator (format " + std::to_string(version) + ")";
        return false;
    }

    // Everything is read into fresh objects first, so a bad file changes nothing.
    Settings s;
    s.vsync = settings.vsync;
    Camera c;
    auto b = std::make_unique<BodySystem>();
    struct Part {
        const char* key;
        std::function<bool(const Value&)> read;
    };
    const Part parts[] = {
        {"settings", [&](const Value& x) { return decodeFields(x, s); }},
        {"camera", [&](const Value& x) { return decodeFields(x, c); }},
        {"bodies", [&](const Value& x) { return decodeFields(x, *b); }},
    };
    for (const Part& part : parts) {
        const Value* x = root.find(part.key);
        g_where.clear();
        if (!x || !part.read(*x)) {
            if (!x) error = std::string("no ") + part.key + " in the file";
            else error = std::string("bad value at ") + part.key + (g_where.empty() ? "" : (g_where[0] == '[' ? "" : ".") + g_where);
            return false;
        }
    }
    if (!checkSlots(*b, error)) return false;

    // The history's samples.
    const std::vector<char>* history = blob.find("HIST");
    std::vector<bool> filled(History::kSlots, false);
    if (history) {
        size_t offset = 0;
        while (offset < history->size()) {
            int32_t slot = -1;
            if (!take(*history, offset, slot) || slot < 0 || slot >= History::kSlots ||
                history->size() - offset < kHistorySlotBytes) {
                error = "the bodies' history in the sidecar is damaged";
                return false;
            }
            std::memcpy(&b->history.image[static_cast<size_t>(History::kRows * slot) * History::kColumns],
                        history->data() + offset, kHistorySlotBytes);
            offset += kHistorySlotBytes;
            filled[static_cast<size_t>(slot)] = true;
        }
    }
    for (int slot = 0; slot < History::kSlots; ++slot) {
        if (b->history.length[slot] > 0 && !filled[static_cast<size_t>(slot)]) {
            error = "the sidecar lacks the history of slot " + std::to_string(slot);
            return false;
        }
    }
    // The renderer re-uploads the whole history when its generation changes.
    b->history.generation = bodies.history.generation + 1;
    b->history.dirtyRows.clear();

    settings = std::move(s);
    camera = c;
    bodies = std::move(*b);
    return true;
}

size_t Blob::bytes() const {
    size_t total = 0;
    for (const auto& chunk : m_chunks) total += chunk.second.size();
    return total;
}

std::string simulationSettings(const Settings& settings) {
    Writer w;
    simulationFields(w, const_cast<Settings&>(settings));
    return json::dump(w.object);
}

bool applySimulationSettings(const std::string& text, Settings& settings) {
    Value root;
    std::string error;
    if (!json::parse(text, root, error) || !root.isObject()) return false;
    Settings s = settings;
    Reader r(root);
    simulationFields(r, s);
    if (!r.ok()) return false;
    settings = std::move(s);
    return true;
}

std::filesystem::path sidecarPath(const std::filesystem::path& path) {
    std::filesystem::path p = path;
    p += ".bin";
    return p;
}

std::string checksum(const std::vector<char>& bytes) {
    uint64_t hash = 14695981039346656037ull;  // FNV-1a offset basis
    for (const char c : bytes) {
        hash ^= static_cast<unsigned char>(c);
        hash *= 1099511628211ull;  // FNV prime
    }
    char text[17];
    std::snprintf(text, sizeof(text), "%016llx", static_cast<unsigned long long>(hash));
    return text;
}

bool writeFiles(const std::filesystem::path& path, json::Value root, const Blob& blob, std::string& error) {
    const std::vector<char> binary = blob.serialize();
    const std::filesystem::path side = sidecarPath(path);
    Value sidecar = Value::object();
    sidecar.set("file", Value::string(utf8(side.filename())));
    sidecar.set("bytes", Value::number(static_cast<long long>(binary.size())));
    sidecar.set("fnv1a", Value::string(checksum(binary)));
    root.set("sidecar", std::move(sidecar));
    const std::string text = json::dump(root);

    std::error_code ec;
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path(), ec);
    std::filesystem::path tmpText = path, tmpSide = side;
    tmpText += ".tmp";
    tmpSide += ".tmp";
    if (!writeWholeFile(tmpSide, binary.data(), binary.size()) || !writeWholeFile(tmpText, text.data(), text.size())) {
        std::filesystem::remove(tmpText, ec);
        std::filesystem::remove(tmpSide, ec);
        error = "could not write " + utf8(path);
        return false;
    }
    std::filesystem::rename(tmpSide, side, ec);
    if (!ec) std::filesystem::rename(tmpText, path, ec);
    if (ec) {
        error = "could not write " + utf8(path) + ": " + ec.message();
        return false;
    }
    return true;
}

bool readFiles(const std::filesystem::path& path, json::Value& root, Blob& blob, std::string& error) {
    std::vector<char> text;
    if (!readWholeFile(path, text)) {
        error = "could not read " + utf8(path);
        return false;
    }
    if (!json::parse(std::string_view(text.data(), text.size()), root, error)) {
        error = utf8(path.filename()) + ": " + error;
        return false;
    }
    blob = Blob();
    const Value* sidecar = root.find("sidecar");
    if (!sidecar) return true;  // written by hand, without one
    const std::filesystem::path side = sidecarPath(path);
    std::vector<char> binary;
    if (!readWholeFile(side, binary)) {
        error = "missing " + utf8(side.filename()) + " (the bodies' past and the gas), which belongs beside the scene";
        return false;
    }
    long long bytes = -1;
    std::string sum;
    const Value* b = sidecar->find("bytes");
    const Value* h = sidecar->find("fnv1a");
    if (!b || !b->get(bytes) || bytes != static_cast<long long>(binary.size()) || !h || !h->get(sum) ||
        sum != checksum(binary)) {
        error = utf8(side.filename()) + " does not belong to this scene (it was changed or saved separately)";
        return false;
    }
    return blob.deserialize(binary, error);
}

}  // namespace scene
