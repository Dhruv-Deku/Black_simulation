#include "keyframes.h"

#include "kerr.h"

#include <algorithm>
#include <cmath>

namespace keyframes {

namespace {

constexpr int kTable = 32;  // arc-length table steps per segment
constexpr double kMinDuration = 0.05;  // seconds: no segment is crossed faster

// Barry and Goldman's pyramid: the Catmull-Rom point at parameter t for points p[0..3]
// with knots k[0..3], t between k[1] and k[2].
glm::dvec3 barryGoldman(const glm::dvec3* p, const double* k, double t) {
    auto lerp = [&](const glm::dvec3& a, const glm::dvec3& b, double ka, double kb) {
        return ((kb - t) * a + (t - ka) * b) / (kb - ka);
    };
    const glm::dvec3 a1 = lerp(p[0], p[1], k[0], k[1]);
    const glm::dvec3 a2 = lerp(p[1], p[2], k[1], k[2]);
    const glm::dvec3 a3 = lerp(p[2], p[3], k[2], k[3]);
    const glm::dvec3 b1 = lerp(a1, a2, k[0], k[2]);
    const glm::dvec3 b2 = lerp(a2, a3, k[1], k[3]);
    return lerp(b1, b2, k[1], k[2]);
}

// Slerp along the arc from a to b as given (squad needs it without the shortcut that
// plain slerp takes when the two lie in opposite hemispheres).
glm::dquat slerpArc(const glm::dquat& a, const glm::dquat& b, double h) {
    const double c = std::clamp(glm::dot(a, b), -1.0, 1.0);
    const double angle = std::acos(c);
    if (angle < 1e-9) return glm::normalize(a * (1.0 - h) + b * h);
    const double s = std::sin(angle);
    return a * (std::sin((1.0 - h) * angle) / s) + b * (std::sin(h * angle) / s);
}

glm::dquat squad(const glm::dquat& q1, const glm::dquat& q2, const glm::dquat& s1, const glm::dquat& s2, double h) {
    return glm::normalize(slerpArc(slerpArc(q1, q2, h), slerpArc(s1, s2, h), 2.0 * h * (1.0 - h)));
}

double smoothstep(double x) {
    x = std::clamp(x, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

// Uniform Catmull-Rom of four values at x in [0, 1] (between b and c).
double catmullRom(double a, double b, double c, double d, double x) {
    return b + 0.5 * x * (c - a + x * (2.0 * a - 5.0 * b + 4.0 * c - d + x * (3.0 * (b - c) + d - a)));
}

glm::dvec3 vec(const std::array<double, 3>& a) {
    return {a[0], a[1], a[2]};
}

// The local r, theta and phi directions as an orthonormal, right-handed Cartesian frame
// (columns). They match the ZAMO's axes up to O(a / r^2), which is how the camera's yaw,
// pitch and roll are read (src/camera.h).
glm::dmat3 localFrame(double spin, double r, double theta, double phi) {
    const kerr::LocalAxes axes = kerr::localAxes(spin, r, theta, phi);
    const glm::dvec3 er = glm::normalize(vec(axes.r));
    glm::dvec3 et = vec(axes.theta);
    et = glm::normalize(et - glm::dot(et, er) * er);
    return glm::dmat3(er, et, glm::cross(er, et));
}

}  // namespace

// --- Quaternions -------------------------------------------------------------------------

glm::dquat logUnit(const glm::dquat& q) {
    const glm::dvec3 v(q.x, q.y, q.z);
    const double n = glm::length(v);
    if (n < 1e-15) return glm::dquat(0.0, 0.0, 0.0, 0.0);
    const double angle = std::atan2(n, q.w);
    const glm::dvec3 axis = v * (angle / n);
    return glm::dquat(0.0, axis.x, axis.y, axis.z);
}

glm::dquat expPure(const glm::dquat& p) {
    const glm::dvec3 v(p.x, p.y, p.z);
    const double angle = glm::length(v);
    if (angle < 1e-15) return glm::dquat(1.0, 0.0, 0.0, 0.0);
    const glm::dvec3 axis = v * (std::sin(angle) / angle);
    return glm::dquat(std::cos(angle), axis.x, axis.y, axis.z);
}

// --- The path ------------------------------------------------------------------------------

void Path::build(const std::vector<Keyframe>& keys) {
    m_keys = keys;
    m_points.clear();
    m_knots.clear();
    m_lengths.clear();
    m_tables.clear();
    m_cumulative.clear();
    m_arrive.clear();
    m_depart.clear();
    m_slope.clear();
    m_quats.clear();
    m_innerOut.clear();
    m_innerIn.clear();
    const size_t n = m_keys.size();
    if (n == 0) return;
    for (Keyframe& k : m_keys) {
        k.duration = std::max(k.duration, kMinDuration);
        k.hold = std::max(k.hold, 0.0);
        k.orientation = glm::normalize(k.orientation);
    }

    // When the camera is at each keyframe.
    m_arrive.resize(n);
    m_depart.resize(n);
    for (size_t i = 0; i < n; ++i) {
        m_arrive[i] = i == 0 ? 0.0 : m_depart[i - 1] + m_keys[i - 1].duration;
        m_depart[i] = m_arrive[i] + m_keys[i].hold;
    }

    // Orientations, each in the hemisphere of the one before (q and -q are the same turn).
    for (size_t i = 0; i < n; ++i) {
        glm::dquat q = m_keys[i].orientation;
        if (i > 0 && glm::dot(q, m_quats.back()) < 0.0) q = -q;
        m_quats.push_back(q);
    }
    m_innerOut = m_innerIn = m_quats;

    m_cumulative.assign(n, 0.0);
    m_slope.assign(n, 0.0);
    if (n < 2) return;

    // The points, with each end reflected through its keyframe, and the centripetal knots:
    // k_j+1 = k_j + |P_j+1 - P_j|^(1/2) (at least a sliver, so that coincident keyframes
    // still have a segment between them).
    m_points.push_back(2.0 * m_keys[0].position - m_keys[1].position);
    for (const Keyframe& k : m_keys) m_points.push_back(k.position);
    m_points.push_back(2.0 * m_keys[n - 1].position - m_keys[n - 2].position);
    double extent = 1.0;
    for (const glm::dvec3& p : m_points) extent = std::max(extent, glm::length(p));
    const double sliver = 1e-6 * std::sqrt(extent);
    m_knots.push_back(0.0);
    for (size_t j = 1; j < m_points.size(); ++j) {
        m_knots.push_back(m_knots.back() + std::max(std::sqrt(glm::length(m_points[j] - m_points[j - 1])), sliver));
    }

    // Arc length along each segment, tabulated.
    for (size_t i = 0; i + 1 < n; ++i) {
        std::vector<double> table(kTable + 1, 0.0);
        for (int k = 0; k < kTable; ++k) {
            table[k + 1] = table[k] + arcLength(i, static_cast<double>(k) / kTable, static_cast<double>(k + 1) / kTable);
        }
        m_lengths.push_back(table.back());
        m_tables.push_back(std::move(table));
        m_cumulative[i + 1] = m_cumulative[i] + m_lengths.back();
    }

    // The speed passing each keyframe: zero where the camera starts, holds or ends;
    // elsewhere the mean of the speeds on either side, limited so that the distance never
    // runs backwards (Fritsch and Carlson: (m_i / d)^2 + (m_i+1 / d)^2 <= 9 on each segment
    // whose mean speed is d).
    std::vector<double> mean(n - 1);
    for (size_t i = 0; i + 1 < n; ++i) mean[i] = m_lengths[i] / (m_arrive[i + 1] - m_depart[i]);
    for (size_t i = 1; i + 1 < n; ++i) {
        if (m_keys[i].hold > 0.0 || mean[i - 1] <= 0.0 || mean[i] <= 0.0) continue;
        m_slope[i] = 0.5 * (mean[i - 1] + mean[i]);
    }
    for (size_t i = 0; i + 1 < n; ++i) {
        if (mean[i] <= 0.0) {
            m_slope[i] = m_slope[i + 1] = 0.0;
            continue;
        }
        const double a = m_slope[i] / mean[i], b = m_slope[i + 1] / mean[i];
        const double r2 = a * a + b * b;
        if (r2 > 9.0) {
            const double tau = 3.0 / std::sqrt(r2);
            m_slope[i] = tau * a * mean[i];
            m_slope[i + 1] = tau * b * mean[i];
        }
    }

    // Squad's inner points. Shoemake's s_i = q_i exp(-(A + B) / 4), with A = log(q_i^-1 q_i+1)
    // and B = log(q_i^-1 q_i-1), turns smoothly only if every segment is crossed at the same
    // rate. Here the orientation follows the distance along the path, so the turning rate
    // per unit length is matched instead: the tangent (A - B) / (L_i-1 + L_i), as in a
    // non-uniform Catmull-Rom, scaled by each segment's length L, gives an inner point for
    // each side, s_out = q_i exp((T_out - A) / 2) and s_in = q_i exp(-(T_in + B) / 2), which
    // reduce to Shoemake's for equal lengths. The ends take their own orientation as the
    // missing neighbour.
    for (size_t i = 0; i < n; ++i) {
        const glm::dquat& q = m_quats[i];
        const glm::dquat inv = glm::conjugate(q);
        const glm::dquat A = i + 1 < n ? logUnit(inv * m_quats[i + 1]) : glm::dquat(0.0, 0.0, 0.0, 0.0);
        const glm::dquat B = i > 0 ? logUnit(inv * m_quats[i - 1]) : glm::dquat(0.0, 0.0, 0.0, 0.0);
        const double before = i > 0 ? m_lengths[i - 1] : 0.0, after = i + 1 < n ? m_lengths[i] : 0.0;
        glm::dquat out = (A - B) * 0.5, in = out;
        if (i > 0 && i + 1 < n && before + after > 0.0) {
            const glm::dquat perLength = (A - B) * (1.0 / (before + after));
            out = perLength * after;
            in = perLength * before;
        }
        m_innerOut[i] = glm::normalize(q * expPure((out - A) * 0.5));
        m_innerIn[i] = glm::normalize(q * expPure((in + B) * -0.5));
    }
}

double Path::duration() const {
    return m_keys.empty() ? 0.0 : m_depart.back();
}

glm::dvec3 Path::curve(size_t segment, double u) const {
    const double* k = &m_knots[segment];
    return barryGoldman(&m_points[segment], k, k[1] + u * (k[2] - k[1]));
}

glm::dvec3 Path::tangent(size_t segment, double u) const {
    const double h = 1e-5;
    return (curve(segment, u + h) - curve(segment, u - h)) / (2.0 * h);
}

double Path::arcLength(size_t segment, double u0, double u1) const {
    // Five-point Gauss-Legendre quadrature of |dP/du|.
    static const double x[5] = {0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640,
                                0.9061798459386640};
    static const double w[5] = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891,
                                0.2369268850561891};
    const double mid = 0.5 * (u0 + u1), half = 0.5 * (u1 - u0);
    double sum = 0.0;
    for (int i = 0; i < 5; ++i) sum += w[i] * glm::length(tangent(segment, mid + half * x[i]));
    return sum * half;
}

double Path::parameterAt(size_t segment, double s) const {
    const std::vector<double>& table = m_tables[segment];
    if (table.back() <= 0.0) return 0.0;
    s = std::clamp(s, 0.0, table.back());
    const int k = std::clamp(static_cast<int>(std::upper_bound(table.begin(), table.end(), s) - table.begin()) - 1, 0,
                             kTable - 1);
    const double u0 = static_cast<double>(k) / kTable, u1 = static_cast<double>(k + 1) / kTable;
    const double span = table[k + 1] - table[k];
    double u = span > 0.0 ? u0 + (u1 - u0) * (s - table[k]) / span : u0;
    // Newton's method on the arc length from the table's step.
    for (int iteration = 0; iteration < 4; ++iteration) {
        const double speed = glm::length(tangent(segment, u));
        if (speed <= 0.0) break;
        u = std::clamp(u - (table[k] + arcLength(segment, u0, u) - s) / speed, u0, u1);
    }
    return u;
}

glm::dvec3 Path::evenly(size_t segment, double s) const {
    return curve(segment, parameterAt(segment, s * m_lengths[segment]));
}

double Path::distance(double seconds, double& rate, size_t& segment) const {
    rate = 0.0;
    const size_t n = m_keys.size();
    for (size_t i = 0; i < n; ++i) {
        segment = std::min(i, n - 2);
        if (seconds <= m_depart[i]) return m_cumulative[i];  // there, or holding there
        if (i + 1 == n) return m_cumulative[i];               // past the end
        if (seconds < m_arrive[i + 1]) {
            // The cubic (Hermite) through both ends at their speeds.
            const double h = m_arrive[i + 1] - m_depart[i];
            const double x = (seconds - m_depart[i]) / h;
            const double x2 = x * x, x3 = x2 * x;
            const double s0 = m_cumulative[i], s1 = m_cumulative[i + 1];
            const double m0 = m_slope[i] * h, m1 = m_slope[i + 1] * h;
            const double s = (2 * x3 - 3 * x2 + 1) * s0 + (x3 - 2 * x2 + x) * m0 + (-2 * x3 + 3 * x2) * s1 + (x3 - x2) * m1;
            rate = ((6 * x2 - 6 * x) * s0 + (3 * x2 - 4 * x + 1) * m0 + (-6 * x2 + 6 * x) * s1 + (3 * x2 - 2 * x) * m1) / h;
            return std::clamp(s, s0, s1);
        }
    }
    return m_cumulative.back();
}

Pose Path::at(double seconds) const {
    Pose pose;
    const size_t n = m_keys.size();
    if (n == 0) return pose;
    if (n == 1) {
        pose.position = m_keys[0].position;
        pose.orientation = m_quats[0];
        pose.fovDeg = m_keys[0].fovDeg;
        return pose;
    }
    double rate = 0.0;
    size_t i = 0;
    const double s = distance(seconds, rate, i);
    // How far along the segment, as a fraction of its length (of its time, if it has no
    // length: the camera then only turns).
    double fraction = 0.0;
    double u = 0.0;
    if (m_lengths[i] > 1e-12 * (1.0 + m_cumulative.back())) {
        u = parameterAt(i, s - m_cumulative[i]);
        fraction = std::clamp((s - m_cumulative[i]) / m_lengths[i], 0.0, 1.0);
        pose.velocity = glm::normalize(tangent(i, u)) * rate;
    } else {
        fraction = smoothstep((seconds - m_depart[i]) / (m_arrive[i + 1] - m_depart[i]));
        u = fraction;
    }
    pose.position = curve(i, u);
    pose.orientation = squad(m_quats[i], m_quats[i + 1], m_innerOut[i], m_innerIn[i + 1], fraction);
    const double f0 = m_keys[i > 0 ? i - 1 : 0].fovDeg, f1 = m_keys[i].fovDeg, f2 = m_keys[i + 1].fovDeg;
    const double f3 = m_keys[std::min(i + 2, n - 1)].fovDeg;
    pose.fovDeg = std::clamp(catmullRom(f0, f1, f2, f3, fraction), 5.0, 150.0);
    return pose;
}

// --- The camera ----------------------------------------------------------------------------

Keyframe fromCamera(const Camera& camera, double spin) {
    Keyframe k;
    const auto p = kerr::toCartesian(spin, camera.r, camera.theta, camera.phi);
    k.position = glm::dvec3(p[0], p[1], p[2]);
    // The camera's right, up and forward axes on the local frame, then in Cartesian. They
    // form a left-handed set (forward is minus right x up), so the rotation is that of
    // right, up and backward.
    const glm::dmat3 local(camera.basis());
    const glm::dmat3 axes = localFrame(spin, camera.r, camera.theta, camera.phi);
    const glm::dmat3 world(axes * local[0], axes * local[1], -(axes * local[2]));
    k.orientation = glm::normalize(glm::quat_cast(world));
    k.fovDeg = camera.fovDeg;
    return k;
}

void applyPose(const Pose& pose, double spin, double minR, Camera& camera) {
    const auto rtp = kerr::fromCartesian(spin, pose.position.x, pose.position.y, pose.position.z);
    camera.r = std::max(rtp[0], minR);
    camera.theta = std::clamp(rtp[1], 0.02, kerr::kPi - 0.02);
    camera.phi = rtp[2];
    const glm::dmat3 world = glm::mat3_cast(glm::normalize(pose.orientation));
    const glm::dmat3 toLocal = glm::transpose(localFrame(spin, camera.r, camera.theta, camera.phi));
    const glm::dvec3 right = toLocal * world[0];
    const glm::dvec3 forward = -(toLocal * world[2]);
    // basis(): forward = (-cos(pitch) cos(yaw), -sin(pitch), cos(pitch) sin(yaw)) on (r, theta, phi);
    // before the roll, right = (sin(yaw), 0, cos(yaw)) and up = right x forward.
    const double pitch = std::asin(std::clamp(-forward.y, -1.0, 1.0));
    const double yaw = std::atan2(forward.z, -forward.x);
    const glm::dvec3 right0(std::sin(yaw), 0.0, std::cos(yaw));
    const glm::dvec3 up0 = glm::cross(right0, glm::normalize(forward));
    camera.pitch = static_cast<float>(pitch);
    camera.yaw = static_cast<float>(yaw);
    camera.roll = static_cast<float>(std::atan2(glm::dot(right, up0), glm::dot(right, right0)));
    camera.fovDeg = static_cast<float>(pose.fovDeg);
}

glm::dvec3 zamoVelocity(const glm::dvec3& position, const glm::dvec3& velocity, double spin, double maxSpeed,
                        bool& capped) {
    capped = false;
    const auto rtp = kerr::fromCartesian(spin, position.x, position.y, position.z);
    const double r = rtp[0], theta = std::clamp(rtp[1], 0.02, kerr::kPi - 0.02);
    if (r < kerr::horizonRadius(spin) + 0.05) {
        capped = true;
        return glm::dvec3(0.0);
    }
    const kerr::Metric m = kerr::metricAt(spin, r, theta);
    // The ZAMO itself is dragged round the axis at omega.
    const glm::dvec3 dragged = m.omega * glm::dvec3(-position.y, position.x, 0.0);
    const double lapse = std::sqrt(m.sigma * m.delta / m.A);
    glm::dvec3 v = glm::transpose(localFrame(spin, r, theta, rtp[2])) * (velocity - dragged) / lapse;
    const double speed = glm::length(v);
    if (speed > maxSpeed) {
        v *= maxSpeed / speed;
        capped = true;
    }
    return v;
}

// --- Files ------------------------------------------------------------------------------------

json::Value toJson(const std::vector<Keyframe>& keys) {
    using json::Value;
    Value list = Value::array();
    for (const Keyframe& k : keys) {
        Value x = Value::object();
        Value p = Value::array(), q = Value::array();
        for (int i = 0; i < 3; ++i) p.push(Value::number(k.position[i]));
        for (double c : {k.orientation.w, k.orientation.x, k.orientation.y, k.orientation.z}) q.push(Value::number(c));
        x.set("position", std::move(p));
        x.set("orientation", std::move(q));  // w, x, y, z
        x.set("fovDeg", Value::number(k.fovDeg));
        x.set("duration", Value::number(k.duration));
        x.set("hold", Value::number(k.hold));
        list.push(std::move(x));
    }
    return list;
}

bool fromJson(const json::Value& value, std::vector<Keyframe>& keys) {
    if (!value.isArray()) return false;
    std::vector<Keyframe> out;
    for (size_t i = 0; i < value.size(); ++i) {
        const json::Value& x = value[i];
        const json::Value *p = x.find("position"), *q = x.find("orientation");
        if (!p || !q || !p->isArray() || !q->isArray() || p->size() != 3 || q->size() != 4) return false;
        Keyframe k;
        double c[4];
        for (int j = 0; j < 3; ++j) {
            if (!(*p)[static_cast<size_t>(j)].get(k.position[j])) return false;
        }
        for (int j = 0; j < 4; ++j) {
            if (!(*q)[static_cast<size_t>(j)].get(c[j])) return false;
        }
        k.orientation = glm::dquat(c[0], c[1], c[2], c[3]);
        const std::pair<const char*, double*> numbers[] = {{"fovDeg", &k.fovDeg}, {"duration", &k.duration}, {"hold", &k.hold}};
        for (const auto& [key, field] : numbers) {
            const json::Value* f = x.find(key);
            if (f && !f->get(*field)) return false;
        }
        out.push_back(k);
    }
    keys = std::move(out);
    return true;
}

}  // namespace keyframes
