#pragma once

// Gravity between orbiting bodies.
//
// Each body follows its Kerr geodesic between kicks; a kick adds the Newtonian pull
// of the other bodies to its 4-velocity, measured in the body's own rest frame:
// each separation runs between events simultaneous for that body, and its length
// comes from the metric. (Measuring them on the coordinate time slices instead, or
// in ZAMO frames, gave fast-moving pairs a spurious steady torque: a moon orbiting a
// planet at 0.2c spun out of its grip within a few orbits.) Kick, drift, kick
// (Strang splitting) makes this second-order accurate in the step. It works
// everywhere, inside the horizon too.

#include "kerr.h"
#include "trajectory.h"

#include <vector>

namespace gravity {

struct Member {
    kerr::Constants c;
    kerr::State s;
    double mass = 0.0;      // in units of the hole's mass
    double softening = 0.1; // closest effective approach (M): bodies have size
};

// The Newtonian pull of the others on member i, as a covariant 4-vector (t, r, theta,
// phi) in i's rest space: m d_mu / |d|^3, with each separation d taken between
// events that are simultaneous for body i and measured with the metric.
kerr::Vec4 pull(double a, const std::vector<Member>& members, size_t i);

// Kicks every member by dt of coordinate time: the pull, projected into the member's
// rest space, changes its 4-velocity over its proper time. This works everywhere,
// inside the horizon too.
void kick(double a, std::vector<Member>& members, double dt);

// Relativistic velocity addition: a velocity u measured in a frame that itself moves
// with velocity v (all relative to the same ZAMO).
void addVelocities(const double v[3], const double u[3], double out[3]);

// A step that resolves the tightest interacting pair: a small fraction of its
// mutual orbital period, and never more than maxStep.
double stepLimit(double a, const std::vector<Member>& members, double maxStep);

// True if any two members pull on each other noticeably.
bool interacting(const std::vector<Member>& members);

// Advances all members together by dt (kick, drift, kick in substeps), recording
// each one's path. Members that fall in or escape stop there.
std::vector<trajectory::Path> predict(double a, std::vector<Member> members, double duration, double minSpacing,
                                      size_t maxPoints = 3000);

}  // namespace gravity
