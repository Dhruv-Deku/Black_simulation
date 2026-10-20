#pragma once

#include "central.h"
#include "tde.h"

// Everything the user can tweak from the control panel.
//
// Sizes and times that lay out the scene (disk radii, escape radius, prediction time,
// time scale) are in scene units: sceneScale() M of length (1 for black holes and
// neutron stars, R/6 for bigger objects) and sceneScale()^1.5 M of time, the
// dynamical time at that size. For a black hole they are simply M.
struct Settings {
    // The central object; `spin` below is the Kerr parameter it implies
    // (central.metricSpin()), kept up to date by the app.
    central::Object central = [] {
        central::Object o;
        o.chi = 0.94;
        return o;
    }();
    bool physicalUnits = false;    // readouts also in km, AU, seconds, hours
    bool radiusFollowsMass = true; // white dwarfs and neutron stars: radius from the mass-radius relation
    bool autoStage = true;         // the compactness slider turns a star into a white dwarf, then a neutron star
    float accretionRate = 0.01f;   // in Eddington units: sets the magnetospheric (Alfven) radius
    bool magnetosphere = true;     // magnetized stars truncate the disk and feed accretion columns
    bool fieldLines = true;        // a neutron star's magnetic field lines, glowing (false color)
    bool pulsarBeams = true;       // pulsars' beams, as if lighting a haze
    bool starquakes = true;        // a young neutron star's crust slips now and then and rings
    float magnetosphereGlow = 1.0f;  // brightness of both
    bool surfaceDetail = true;     // spots, granulation, hot spots, flares
    float surfaceBrightness = 1.0f;
    bool corona = true;            // stellar coronae and planetary atmospheres
    float coronaBoost = 1e5f;      // the corona's brightness, times the real one (a millionth of the disk)
    bool companionStar = true;     // a white dwarf with a disk: its Roche-lobe-filling companion and stream
    bool filmColor = true;
    bool planetMaps = true;        // the Earth and Jupiter from spacecraft maps (assets/planets)         // stars in the warm colors of films and colorized solar images (else true color)

    // Black hole
    float spin = 0.94f;  // a = J/M, in [0, 0.999]

    // Accretion disk
    bool diskEnabled = true;
    bool diskPrograde = true;      // gas orbits with the hole's rotation
    bool diskInnerAtIsco = true;   // inner edge follows the ISCO as spin changes
    float diskInner = 2.0f;        // used when diskInnerAtIsco is off
    float diskOuter = 20.0f;
    float diskTemperature = 6500.0f;  // peak temperature in K
    // Around a neutron star or white dwarf the app sets the disk's real peak temperature
    // from the accretion rate (K; 0 otherwise), and diskTemperature to its display value.
    double diskPhysicalTemperature = 0.0;
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
    // Display only: the debris glows at least as if it were this hot (0: as hot as it is).
    // A disrupted star's stream is lit by its own leftover heat and then by the flare.
    float debrisGlowFloor = 0.0f;

    // After a tidal disruption (src/tde.h): the debris that comes back powers a flare,
    // feeds a small new disk and blows out a billowing outflow along its orbit's axis.
    bool tdeFlare = true;
    float outflowFraction = 0.5f;   // share of the returning debris blown out (the rest is accreted)
    float outflowSpeed = 0.1f;      // terminal speed on the axis, in c (AT2019qiz: ~0.03)
    float outflowOpeningDeg = 14.0f;  // half-opening angle once the fallback has faded
    float outflowWidening = 3.0f;   // ... times this at first, before the new disk collimates it
    float outflowThickness = 1.0f;  // opacity (a real wind is far more opaque; see README)
    float outflowBrightness = 1.0f; // how brightly the flare lights it
    bool flareDisk = true;          // the returning gas forms a small disk that shines with the flare
    bool flareCinematic = true;     // the flare light fades as R^-0.7, not R^-2, so the outflow's top shows
    bool debrisFilaments = false;   // the star's structure, stretched by the tide into fine strands (costly)
    bool flareFastForward = true;   // time runs faster while the bound debris is far out, until it returns
    bool flareCamera = false;       // the camera's distance follows a tidal disruption (until you zoom)

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
    // Traced below the window's resolution, the image is enlarged with FSR 1 (an
    // edge-adaptive upscale, then contrast-adaptive sharpening) instead of bilinearly.
    int upscaler = 1;            // 0 = bilinear, 1 = FSR 1
    float sharpness = 0.25f;     // RCAS, in stops below the sharpest
    // The resolution follows a frame-time target: scale *= (target / measured)^(1/2), smoothed.
    bool dynamicResolution = false;
    float targetMs = 16.7f;

    // Simulation
    float timeScale = 15.0f;  // M of coordinate time per real second
    bool penroseSpinDown = true;  // a Penrose split's negative-energy fragment spins the hole down
    bool paused = false;

    // Sound (src/scene_sound.h): off by default, and quiet when on.
    bool sound = false;
    float soundVolume = 0.5f;
    bool soundBodies = true;    // a tone for each body
    bool soundPulsars = true;   // pulsars' clicks
    bool soundRumble = true;    // the roar of the tide near the hole
    bool soundWaves = true;     // gravitational waves of heavy bodies, raised in pitch

    // Tools
    bool rayInspector = true;  // Ctrl+click or middle-click traces the light reaching a pixel

    // Camera
    bool slowNearHole = true;  // slow time down while the camera falls in close to the hole
    bool autoOrbit = false;
    float autoOrbitSpeed = 0.06f;  // rad/s
};

// How a tidal disruption's outflow is launched, from the panel's settings.
inline tde::OutflowParams outflowParams(const Settings& s) {
    tde::OutflowParams p;
    p.fraction = s.outflowFraction;
    p.speed = s.outflowSpeed;
    p.openingLate = s.outflowOpeningDeg * 3.14159265358979 / 180.0;
    p.widening = s.outflowWidening;
    return p;
}
