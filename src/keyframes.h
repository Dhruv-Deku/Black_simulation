#pragma once

// Cinematic camera paths through keyframes.
//
// Position: a centripetal Catmull-Rom spline (alpha = 1/2) through the keyframes'
// positions, evaluated with Barry and Goldman's pyramid of linear interpolations. It
// passes through every keyframe and, unlike the uniform kind, never forms a cusp or a
// loop inside a segment (Yuksel, Schaefer and Keyser 2011). Each end is extended by
// reflecting its neighbour.
// Speed: each segment is reparametrized by arc length (Gauss-Legendre quadrature, then
// inverted by Newton's method), so the camera moves evenly along it. How far along the path
// it is follows a monotone cubic (Fritsch and Carlson 1980) through the times it is due
// at each keyframe: at rest where it starts, holds and ends, smooth in between.
// Orientation: squad (Shoemake 1987) between the keyframes' quaternions, its inner points
// weighted by the segments' lengths so that it turns smoothly through them. Field of
// view: a Catmull-Rom of the keyframes' values.
//
// Positions are Kerr-Schild Cartesian, in M; orientations rotate the camera's own axes
// (right, up, forward) onto Cartesian ones.

#include "camera.h"
#include "json.h"

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include <string>
#include <vector>

namespace keyframes {

struct Keyframe {
    glm::dvec3 position{0.0};
    glm::dquat orientation{1.0, 0.0, 0.0, 0.0};
    double fovDeg = 60.0;
    double duration = 4.0;  // seconds from here to the next keyframe
    double hold = 0.0;      // seconds the camera stays here first
};

struct Pose {
    glm::dvec3 position{0.0};
    glm::dquat orientation{1.0, 0.0, 0.0, 0.0};
    double fovDeg = 60.0;
    glm::dvec3 velocity{0.0};  // M per second along the path
};

class Path {
public:
    void build(const std::vector<Keyframe>& keys);
    bool empty() const { return m_keys.empty(); }
    // Seconds from the start to the end of the last keyframe's hold.
    double duration() const;
    Pose at(double seconds) const;

    size_t segments() const { return m_keys.size() < 2 ? 0 : m_keys.size() - 1; }
    // Segment i (keyframe i to i + 1) at spline parameter u in [0, 1] ...
    glm::dvec3 curve(size_t segment, double u) const;
    // ... and at fraction s of its arc length.
    glm::dvec3 evenly(size_t segment, double s) const;
    double segmentLength(size_t segment) const { return m_lengths[segment]; }
    // When the camera arrives at keyframe i.
    double arrival(size_t i) const { return m_arrive[i]; }

private:
    double arcLength(size_t segment, double u0, double u1) const;
    double parameterAt(size_t segment, double s) const;  // u at arc length s (not a fraction)
    glm::dvec3 tangent(size_t segment, double u) const;  // dP/du
    // Distance along the path at a time, and its rate; which segment that is in.
    double distance(double seconds, double& rate, size_t& segment) const;

    std::vector<Keyframe> m_keys;
    std::vector<glm::dvec3> m_points;  // the keyframes' positions, with one more at each end
    std::vector<double> m_knots;       // their centripetal parameters
    std::vector<double> m_lengths;     // of each segment
    std::vector<std::vector<double>> m_tables;  // arc length at kTable + 1 even steps of u, per segment
    std::vector<double> m_cumulative;  // arc length from the start to each keyframe
    std::vector<double> m_arrive, m_depart;  // times at each keyframe (before and after its hold)
    std::vector<double> m_slope;       // the path's speed (M per second) passing each keyframe
    std::vector<glm::dquat> m_quats;  // orientations, each in the hemisphere of the last
    std::vector<glm::dquat> m_innerOut, m_innerIn;  // squad's inner points leaving and reaching each keyframe
};

// The keyframes being edited, and the path through them.
struct Track {
    std::vector<Keyframe> keys;
    Path path;
    bool playing = false;
    double time = 0.0;      // seconds along the path
    bool physical = false;  // move for real (relativistic aberration), not hover at each point
    bool aim = false;       // keep facing the center between keyframes (the keyframes' roll is kept)
    std::string status;     // what playing it last ran into
    void rebuild() { path.build(keys); }
};

// The camera as a keyframe (it must be hovering or moving outside the horizon), and a
// pose put onto a hovering camera: position (kept above minR) and orientation, as yaw,
// pitch and roll in its local frame (src/camera.h).
Keyframe fromCamera(const Camera& camera, double spin);
void applyPose(const Pose& pose, double spin, double minR, Camera& camera);

// The camera's 3-velocity relative to the local zero-angular-momentum observer for a path
// moving at velocity (Kerr-Schild Cartesian, M per M of coordinate time) through position:
// the coordinate velocity on the local axes, less the ZAMO's own frame-dragged motion,
// over the lapse (exact far from the hole; near it, the spatial metric's departure from
// flat is left out). Capped at maxSpeed (c = 1); `capped` says if it had to be.
glm::dvec3 zamoVelocity(const glm::dvec3& position, const glm::dvec3& velocity, double spin, double maxSpeed,
                        bool& capped);

// For scene files.
json::Value toJson(const std::vector<Keyframe>& keys);
bool fromJson(const json::Value& value, std::vector<Keyframe>& keys);

// Quaternion log and exp (unit and pure quaternions), as squad uses them.
glm::dquat logUnit(const glm::dquat& q);
glm::dquat expPure(const glm::dquat& p);

}  // namespace keyframes
