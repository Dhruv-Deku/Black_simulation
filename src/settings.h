#pragma once

// Everything the user can tweak from the control panel.
struct Settings {
    // Black hole
    float spin = 0.94f;  // a = J/M, in [0, 0.999]

    // Accretion disk
    bool diskEnabled = true;
    bool diskPrograde = true;      // gas orbits with the hole's rotation
    bool diskInnerAtIsco = true;   // inner edge follows the ISCO as spin changes
    float diskInner = 2.0f;        // used when diskInnerAtIsco is off
    float diskOuter = 20.0f;
    float diskTemperature = 6500.0f;  // peak temperature in K
    float diskBrightness = 4.0f;
    float diskOpacity = 0.95f;
    float diskTurbulence = 0.75f;
    int relativity = 0;  // 0 = full (Doppler + gravitational shift), 1 = beaming only, 2 = off
    bool lightDelay = true;  // draw things where they were when their light left them
    bool bodyGravity = true; // bodies pull on each other (and on debris gas)

    // Debris gas (tidal disruptions, clouds, rings)
    float debrisBrightness = 2.5f;  // emission per unit density
    float debrisOpacity = 0.1f;     // absorption per unit density
    bool gasPhysics = true;         // pressure, shocks, heating and cooling
    float gasCollisions = 0.3f;     // how readily crossing streams shock (per M per unit density)
    float gasCoolingTime = 60.0f;   // M, at 10,000 K
    float gasFloorTemperature = 4500.0f;

    // Relativistic jets along the spin axis
    bool jetsEnabled = false;
    float jetPower = 1.0f;
    float jetSpeed = 0.8f;       // v / c
    float jetOpeningDeg = 5.0f;  // half opening angle

    // Orbit paths drawn into the scene (ray traced, so they are lensed too)
    bool showPaths = true;
    bool showPredictions = true;
    float predictionTime = 400.0f;  // M of coordinate time ahead
    float pathWidth = 1.3f;         // pixels

    // Sky
    int skyMode = 0;  // 0 = stars, 1 = grid, 2 = image, 3 = none
    float starBrightness = 1.0f;
    float galaxyBrightness = 1.0f;
    float skyImageBrightness = 1.0f;

    // Rendering
    float renderScale = 0.75f;
    int maxSteps = 700;
    float stepSize = 0.04f;
    float escapeRadius = 250.0f;
    float exposure = 1.0f;
    float bloomStrength = 0.05f;
    int tonemapper = 0;  // 0 = ACES, 1 = Reinhard, 2 = linear clamp
    int debugView = 0;   // 0 = off, 1 = step count heat map
    bool vsync = true;

    // Simulation
    float timeScale = 15.0f;  // M of coordinate time per real second
    bool paused = false;

    // Camera
    bool slowNearHole = true;  // slow time down while the camera falls in close to the hole
    bool autoOrbit = false;
    float autoOrbitSpeed = 0.06f;  // rad/s
};
