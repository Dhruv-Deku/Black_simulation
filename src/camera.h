#pragma once

#include "kerr.h"

#include <glm/glm.hpp>

#include <array>
#include <string>

enum class CameraMotion {
    Hover,      // held in place by a zero-angular-momentum observer (ZAMO)
    Orbit,      // riding a near-circular orbit
    FreeFall,   // dropped from rest
    RideAlong,  // sitting on one of the orbiting bodies
};

// A camera at position (r, theta, phi). It looks at the hole by default; yaw and
// pitch turn the view away from that direction. When it moves, it follows its own
// geodesic and sees the sky through relativistic aberration. A falling camera
// crosses the horizon and keeps going until just above the inner horizon; above a
// star or planet it lands on the surface and stays there.
struct Camera {
    double r = 32.0;
    double theta = 1.48;  // polar angle from the spin axis
    double phi = 0.0;
    float yaw = 0.0f;
    float pitch = 0.0f;
    float roll = 0.0f;  // turns the picture about the view direction (radians)
    float fovDeg = 60.0f;

    CameraMotion motion = CameraMotion::Hover;
    kerr::Constants geodesic;  // the path while moving
    kerr::State path;
    double properTime = 0.0;   // time elapsed on the camera's own clock while moving
    bool ended = false;        // the fall reached the interior limit; the camera stops there
    int rideSlot = -1;         // history slot of the body being ridden
    bool lookAhead = false;    // while moving, face the direction of travel

    // The central object, as the camera needs it (kept up to date by the app).
    double surface = 0.0;      // surface radius in M; 0 for a black hole
    double scale = 1.0;        // scene scale: how far out the camera may go
    bool landed = false;       // it came down on the surface (cleared by the app once reported)
    double landingSpeed = 0.0; // how fast it was moving, relative to a local static observer
    // A hovering camera carried along a keyframed path for real (src/keyframes.h): its
    // velocity relative to the local ZAMO, on the r, theta and phi axes (zero otherwise).
    // It sees the sky aberrated and Doppler shifted, as a moving camera does.
    std::array<double, 3> drift{};

    // Columns are the camera's right, up and forward axes, written on the frame's
    // spatial axes (which follow r, theta and phi).
    glm::mat3 basis() const;
    // The camera's orthonormal frame: its 4-velocity and spatial axes.
    kerr::Frame frame(double spin) const;
    bool insideHorizon(double spin) const { return surface <= 0.0 && r < kerr::horizonRadius(spin); }
    // The closest the camera may hover: just above the horizon, or the surface.
    double minDistance(double spin) const;

    void orbit(double dPhi, double dTheta);
    void zoom(double factor, double minR);
    void look(float dYaw, float dPitch);
    void reset();

    // Start or stop moving. Returns false with a reason if that is impossible here.
    bool startOrbit(double spin, std::string& message);
    bool startFall(double spin, std::string& message);
    bool hover(double spin, std::string& message);
    // Leaves the camera on a geodesic (when the body it rode is gone, for example).
    void coast(const kerr::Constants& c, const kerr::State& s);
    // Moves the camera onto a body's path (called every frame while riding).
    void rideWith(int slot, const kerr::Constants& c, const kerr::State& s, double spin, double dt);
    // Moves along the path by dt of coordinate time.
    void advance(double spin, double dt);
    void respin(double oldSpin, double newSpin);

    // Velocity relative to the local ZAMO (zero while hovering; zero inside the
    // horizon, where there are no ZAMOs).
    glm::vec3 velocity(double spin) const;
    // Turns the view to face the direction of travel.
    void faceVelocity(double spin);
    // How fast the camera's clock runs compared with coordinate time (dtau/dt).
    double clockRate(double spin) const;
};
