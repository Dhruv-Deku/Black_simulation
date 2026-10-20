#include "camera.h"

#include <algorithm>
#include <cmath>

glm::mat3 Camera::basis() const {
    // Facing the hole: forward = -e_r, up = -e_theta (north), right = e_phi.
    const glm::vec3 forward0(-1.0f, 0.0f, 0.0f);
    const glm::vec3 up0(0.0f, -1.0f, 0.0f);
    const glm::vec3 right0(0.0f, 0.0f, 1.0f);

    const float cy = std::cos(yaw), sy = std::sin(yaw);
    const float cp = std::cos(pitch), sp = std::sin(pitch);
    const glm::vec3 yawed = cy * forward0 + sy * right0;
    const glm::vec3 right = cy * right0 - sy * forward0;
    const glm::vec3 forward = cp * yawed + sp * up0;
    const glm::vec3 up = glm::cross(right, forward);
    // Roll turns right and up about the view direction.
    const float cr = std::cos(roll), sr = std::sin(roll);
    return glm::mat3(cr * right + sr * up, cr * up - sr * right, forward);
}

kerr::Frame Camera::frame(double spin) const {
    if (motion == CameraMotion::Hover) {
        const kerr::Vec4 zamo = kerr::zamoFourVelocity(spin, r, theta);
        const double v2 = drift[0] * drift[0] + drift[1] * drift[1] + drift[2] * drift[2];
        if (v2 <= 0.0 || v2 >= 1.0 || insideHorizon(spin)) return kerr::observerFrame(spin, r, theta, zamo);
        // Boosted from the ZAMO: u = gamma (e0 + v^i e_i).
        const kerr::Frame z = kerr::observerFrame(spin, r, theta, zamo);
        const double gamma = 1.0 / std::sqrt(1.0 - v2);
        kerr::Vec4 u;
        for (int i = 0; i < 4; ++i) {
            u[i] = gamma * (z.e[0][i] + drift[0] * z.e[1][i] + drift[1] * z.e[2][i] + drift[2] * z.e[3][i]);
        }
        return kerr::observerFrame(spin, r, theta, u);
    }
    kerr::State s = path;
    s.theta = theta;  // clamped off the axis, like the position
    return kerr::observerFrame(spin, r, theta, kerr::fourVelocity(spin, geodesic, s));
}

void Camera::orbit(double dPhi, double dTheta) {
    phi = std::fmod(phi + dPhi, 2.0 * kerr::kPi);
    theta = std::clamp(theta + dTheta, 0.02, kerr::kPi - 0.02);
}

void Camera::zoom(double factor, double minR) {
    r = std::clamp(r * factor, minR, 2000.0 * scale);
}

double Camera::minDistance(double spin) const {
    return surface > 0.0 ? surface + 0.02 * scale : kerr::horizonRadius(spin) + 0.3;
}

void Camera::look(float dYaw, float dPitch) {
    yaw = std::remainder(yaw + dYaw, 6.2831853f);
    pitch = std::clamp(pitch + dPitch, -1.5f, 1.5f);
}

void Camera::reset() {
    *this = Camera();
}

bool Camera::startOrbit(double spin, std::string& message) {
    if (insideHorizon(spin)) {
        message = "Nothing can orbit inside the horizon.";
        return false;
    }
    if (r < kerr::photonOrbitRadius(spin, true) * 1.05) {
        message = "Too close to orbit: inside the photon orbit nothing can circle the hole.";
        return false;
    }
    // Along e_phi at the current latitude: an orbit tilted by (90 deg - theta).
    const double v[3] = {0.0, 0.0, kerr::circularZamoSpeed(spin, r, true)};
    kerr::launchFromZamo(spin, r, theta, phi, v, 1.0, geodesic, path);
    motion = CameraMotion::Orbit;
    properTime = 0.0;
    ended = false;
    message = r < kerr::iscoRadius(spin, true) ? "Orbiting inside the ISCO: the orbit is unstable."
                                                : "Orbiting.";
    return true;
}

bool Camera::startFall(double spin, std::string& message) {
    if (motion != CameraMotion::Hover) {
        // Already moving: keep the current path, just stop following a body.
        motion = CameraMotion::FreeFall;
        rideSlot = -1;
        message = "Falling freely.";
        return true;
    }
    if (surface <= 0.0 && r < kerr::horizonRadius(spin) + 0.2) {
        message = "Too close to the horizon to start from rest.";
        return false;
    }
    if (surface > 0.0 && r < surface + 0.05 * scale) {
        message = "Already on the surface.";
        return false;
    }
    const double v[3] = {0.0, 0.0, 0.0};  // at rest relative to the ZAMO: zero angular momentum
    kerr::launchFromZamo(spin, r, theta, phi, v, 1.0, geodesic, path);
    motion = CameraMotion::FreeFall;
    properTime = 0.0;
    ended = false;
    message = "Falling.";
    return true;
}

bool Camera::hover(double spin, std::string& message) {
    if (surface <= 0.0 && r < kerr::horizonRadius(spin) + 0.05) {
        message = "Nothing can hover inside the horizon: everything there must fall inward.";
        return false;
    }
    motion = CameraMotion::Hover;
    ended = false;
    rideSlot = -1;
    message = "Hovering.";
    return true;
}

void Camera::coast(const kerr::Constants& c, const kerr::State& s) {
    motion = CameraMotion::FreeFall;
    rideSlot = -1;
    geodesic = c;
    path = s;
}

void Camera::rideWith(int slot, const kerr::Constants& c, const kerr::State& s, double spin, double dt) {
    if (motion != CameraMotion::RideAlong || rideSlot != slot) {
        motion = CameraMotion::RideAlong;
        rideSlot = slot;
        properTime = 0.0;
        ended = false;
    } else if (dt > 0.0) {
        properTime += dt * clockRate(spin);
    }
    geodesic = c;
    path = s;
    r = s.r;
    theta = std::clamp(s.theta, 0.02, kerr::kPi - 0.02);
    phi = std::fmod(s.phi, 2.0 * kerr::kPi);
}

void Camera::advance(double spin, double dt) {
    if (motion == CameraMotion::Hover || motion == CameraMotion::RideAlong || ended || dt <= 0.0) return;
    const double limit = kerr::interiorLimit(spin);
    const double ground = surface + 0.02 * scale;  // where a landing camera stops (stars and planets)
    const double target = path.t + dt;
    for (int guard = 0; path.t < target - 1e-12 && guard < 20000; ++guard) {
        const double h = std::min(target - path.t, kerr::suggestedTimeStep(spin, geodesic, path, 0.005));
        const kerr::Derivatives d = kerr::derivatives(spin, geodesic, path);
        properTime += h * kerr::metricAt(spin, path.r, path.theta).sigma / d.dt;  // dtau/dt = Sigma / (dt/dMino)
        const kerr::State before = path;
        kerr::stepCoordinateTime(spin, geodesic, path, h);
        kerr::normalizeAngles(path);
        if (surface > 0.0 && path.r <= ground) {
            // Touchdown: measure the speed against the local ZAMO just before, then stand
            // on the surface (held there like a hovering observer).
            double v[3];
            kerr::zamoVelocity(spin, geodesic, before, v);
            landingSpeed = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            landed = true;
            motion = CameraMotion::Hover;
            rideSlot = -1;
            r = ground;
            theta = std::clamp(before.theta, 0.02, kerr::kPi - 0.02);
            phi = std::fmod(before.phi, 2.0 * kerr::kPi);
            // Look out along the ground toward the horizon (on a compact star it is lifted
            // by light bending).
            yaw = 0.0f;
            pitch = 1.4f;  // nearly level, the ground just below the view
            lookAhead = false;
            return;
        }
        if ((surface <= 0.0 && path.r < limit) || !std::isfinite(path.r)) {
            ended = true;
            break;
        }
    }
    r = path.r;
    theta = std::clamp(path.theta, 0.02, kerr::kPi - 0.02);
    phi = std::fmod(path.phi, 2.0 * kerr::kPi);
}

void Camera::respin(double oldSpin, double newSpin) {
    if (motion == CameraMotion::Hover || motion == CameraMotion::RideAlong) return;
    const double t = path.t;
    if (path.r > kerr::horizonRadius(oldSpin) + 0.05 && path.r > kerr::horizonRadius(newSpin) + 0.05) {
        double v[3];
        kerr::zamoVelocity(oldSpin, geodesic, path, v);
        kerr::launchFromZamo(newSpin, path.r, path.theta, path.phi, v, 1.0, geodesic, path);
    } else {
        // No ZAMO to measure against: keep the constants and re-solve the velocities.
        path.vr = -std::sqrt(std::max(kerr::radialPotential(newSpin, geodesic, path.r), 0.0));
        path.vtheta = std::copysign(std::sqrt(std::max(kerr::polarPotential(newSpin, geodesic, path.theta), 0.0)),
                                    path.vtheta);
    }
    path.t = t;
}

glm::vec3 Camera::velocity(double spin) const {
    if (motion == CameraMotion::Hover) {
        return glm::vec3(static_cast<float>(drift[0]), static_cast<float>(drift[1]), static_cast<float>(drift[2]));
    }
    if (surface <= 0.0 && r < kerr::horizonRadius(spin) + 0.05) return glm::vec3(0.0f);
    double v[3];
    kerr::zamoVelocity(spin, geodesic, path, v);
    return glm::vec3(static_cast<float>(v[0]), static_cast<float>(v[1]), static_cast<float>(v[2]));
}

void Camera::faceVelocity(double spin) {
    const glm::vec3 v = velocity(spin);
    if (glm::length(v) < 1e-4f) return;
    // In the camera's frame its direction of travel has the same components as its
    // velocity relative to the ZAMO (for a boost along r or phi exactly).
    // basis(): forward = (-cos(pitch) cos(yaw), -sin(pitch), cos(pitch) sin(yaw)) in (r, theta, phi).
    const glm::vec3 d = glm::normalize(v);
    pitch = std::asin(std::clamp(-d.y, -1.0f, 1.0f));
    yaw = std::atan2(d.z, -d.x);
}

double Camera::clockRate(double spin) const {
    if (motion == CameraMotion::Hover) {
        // A ZAMO's clock runs at sqrt(Sigma Delta / A) of a distant one.
        const kerr::Metric m = kerr::metricAt(spin, r, theta);
        return std::sqrt(m.sigma * m.delta / m.A);
    }
    const kerr::Derivatives d = kerr::derivatives(spin, geodesic, path);
    return kerr::metricAt(spin, path.r, path.theta).sigma / d.dt;
}
