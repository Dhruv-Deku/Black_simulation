// The surface of a star or planet at the center, its air or corona.
// (Included by trace.frag after the Kerr helpers and shaders/surfaces.glsl.)

const int kLookNone = 0, kLookPhotosphere = 1, kLookBrownDwarf = 2, kLookWhiteDwarf = 3, kLookNeutronStar = 4,
          kLookEarth = 5, kLookGasGiant = 6;
// Each look is compiled into its own variant of the shader (src/renderer.cpp): the
// registers an Earth's sky and continents need would otherwise slow a neutron star's
// ray tracing too.
#ifndef SURFACE_LOOK
#define SURFACE_LOOK 0
#endif

// Brightness of a star's disk center (before exposure), chosen so the tone curve keeps
// its color, limb darkening and spots instead of flattening it to white.
const float kPhotosphereLevel = 0.42;

// --- Film colors -------------------------------------------------------------------------
//
// A blackbody's true color is pale: the Sun is white, a 3,000 K red dwarf peach. Films and
// colorized solar images (SDO, SOHO) show stars far warmer and more saturated: the Sun
// golden with an orange limb, red dwarfs red-orange, brown dwarfs deep red. With
// uFilmColor only the hue is changed: brightness keeps the physics (limb darkening, spots,
// granulation, flares, Doppler and gravitational shifts). (uFilmColor: trace.frag.)
// Brightness of the disk center in the film palette (it peaks toward white, as an
// exposure set for the limb leaves the center).
const float kFilmLevel = 0.95;
// Brightness of the 30.4 nm color table's white.
const float kSdoLevel = 1.15;

// The film hue of temperature T, its brightest channel 1: stops in log10 T.
vec3 filmHue(float T) {
    const int N = 8;
    const float stops[N] = float[N](2.9, 3.18, 3.48, 3.6, 3.76, 3.88, 4.0, 4.5);
    const vec3 hues[N] = vec3[N](vec3(0.5, 0.04, 0.06), vec3(0.92, 0.12, 0.03), vec3(1.0, 0.17, 0.02),
                                 vec3(1.0, 0.33, 0.04), vec3(1.0, 0.56, 0.12), vec3(1.0, 0.84, 0.55),
                                 vec3(0.86, 0.9, 1.0), vec3(0.62, 0.74, 1.0));
    float l = log2(max(T, 1.0)) * 0.30103;
    vec3 c = hues[N - 1];
    if (l <= stops[0]) {
        c = hues[0];
    } else {
        for (int k = 1; k < N; ++k) {
            if (l <= stops[k]) {
                c = mix(hues[k - 1], hues[k], (l - stops[k - 1]) / (stops[k] - stops[k - 1]));
                break;
            }
        }
    }
    return c / max(c.r, max(c.g, c.b));
}

// Thermal light in the film palette: the shifted blackbody's luminance (relative to Tref)
// in the hue of the temperature seen.
vec3 filmBlackbody(float T, float g, float Tref) {
    return filmHue(T * g) * (luminance(shiftedBlackbody(T, g, Tref)) * kFilmLevel);
}

// The electric blues films and artists' impressions give neutron stars, magnetars and
// white dwarfs (their light is mostly ultraviolet and X-rays): deep blue, blue, cyan,
// white; past 1 brighter white still (it blooms).
vec3 electricHue(float t) {
    const int N = 5;
    const float stops[N] = float[N](0.0, 0.3, 0.55, 0.8, 1.0);
    const vec3 cols[N] = vec3[N](vec3(0.0, 0.01, 0.06), vec3(0.03, 0.12, 0.6), vec3(0.1, 0.45, 1.0),
                                 vec3(0.55, 0.88, 1.0), vec3(1.0));
    t = max(t, 0.0);
    if (t >= 1.0) return vec3(1.0 + 1.5 * (t - 1.0));
    for (int k = 1; k < N; ++k) {
        if (t <= stops[k]) return mix(cols[k - 1], cols[k], (t - stops[k - 1]) / (stops[k] - stops[k - 1]));
    }
    return cols[N - 1];
}
const float kElectricLevel = 1.5;

// The saturated blues of films' spinning neutron stars: deep blue, blue, cyan, then
// white only at the very top (a glowing rim).
vec3 beaconHue(float t) {
    const int N = 5;
    const float stops[N] = float[N](0.0, 0.35, 0.65, 0.85, 1.0);
    const vec3 cols[N] = vec3[N](vec3(0.0, 0.03, 0.22), vec3(0.0, 0.28, 0.9), vec3(0.01, 0.68, 1.0),
                                 vec3(0.35, 0.92, 1.0), vec3(1.0));
    t = max(t, 0.0);
    if (t >= 1.0) return vec3(1.0 + 1.5 * (t - 1.0));
    for (int k = 1; k < N; ++k) {
        if (t <= stops[k]) return mix(cols[k - 1], cols[k], (t - stops[k - 1]) / (stops[k] - stops[k - 1]));
    }
    return cols[N - 1];
}

// SDO's 30.4 nm color table (black, dark red, orange, yellow, white), in which films and
// NASA show the Sun; past white, brighter still (it blooms).
vec3 aia304(float t) {
    const int N = 6;
    const float stops[N] = float[N](0.0, 0.2, 0.42, 0.64, 0.85, 1.1);
    const vec3 cols[N] = vec3[N](vec3(0.0), vec3(0.3, 0.02, 0.0), vec3(0.78, 0.15, 0.0), vec3(1.0, 0.45, 0.04),
                                 vec3(1.0, 0.78, 0.3), vec3(1.0, 0.98, 0.85));
    t = max(t, 0.0);
    if (t >= stops[N - 1]) return cols[N - 1] * (1.0 + 1.5 * (t - stops[N - 1]));
    for (int k = 1; k < N; ++k) {
        if (t <= stops[k]) return mix(cols[k - 1], cols[k], (t - stops[k - 1]) / (stops[k] - stops[k - 1]));
    }
    return cols[N - 1];
}

// The magnetic axis in Kerr-Schild Cartesian space when the star had turned through
// `turn` (its rotation angle then, from phi_BL, plus the Kerr-Schild azimuth offset).
vec3 magneticAxisAt(float turn) {
    float c = cos(turn), s = sin(turn);
    return vec3(c * uMagneticAxis.x - s * uMagneticAxis.y, s * uMagneticAxis.x + c * uMagneticAxis.y, uMagneticAxis.z);
}

#if SURFACE_LOOK == 4
// --- A young neutron star's starquakes ----------------------------------------------------
//
// After a glitch the crust rings in its torsional modes (horizontal shear, at ~22 Hz and
// up; central::crustTorsionalFrequency). The field is anchored in the crust, so its
// footpoints, the magnetic poles, swing back and forth: the beams wobble. For each quake:
// how long ago it began, the poles' swing at first (radians), the mode's angular
// frequency and decay time; which way it shakes them (rotating frame) and its degree l.
uniform int uQuakeCount;
uniform vec4 uQuakes[4];
uniform vec4 uQuakeAxes[4];
uniform int uBeacon;       // a young neutron star drawn as films show one (central::Object::beacon)
uniform float uHaloGlow;   // the glow around it (0: none)

// How far the magnetic poles were displaced (rotating frame) `delay` ago.
vec3 quakeShift(float delay) {
    vec3 shift = vec3(0.0);
    for (int k = 0; k < uQuakeCount; ++k) {
        float age = uQuakes[k].x - delay;
        if (age < 0.0) continue;  // its light has not reached us yet
        shift += uQuakeAxes[k].xyz * (uQuakes[k].y * exp(-age / uQuakes[k].w) * sin(uQuakes[k].z * age));
    }
    return shift;
}

// The magnetic axis (Kerr-Schild Cartesian) when the star had turned through `turn`,
// shaken by the starquakes ringing then (`delay` ago).
vec3 magneticAxisShaken(float turn, float delay) {
    vec3 m = uQuakeCount > 0 ? normalize(uMagneticAxis + quakeShift(delay)) : uMagneticAxis;
    float c = cos(turn), s = sin(turn);
    return vec3(c * m.x - s * m.y, s * m.x + c * m.y, m.z);
}

// The ringing crust is heated where it shears: its surface brightens and dims in the
// mode's pattern, a Legendre polynomial P_l about the direction it shakes.
float quakeShimmer(vec3 n, float delay) {
    float f = 1.0;
    for (int k = 0; k < uQuakeCount; ++k) {
        float age = uQuakes[k].x - delay;
        if (age < 0.0) continue;
        float swing = uQuakes[k].y * exp(-age / uQuakes[k].w) * sin(uQuakes[k].z * age);
        float x = dot(n, uQuakeAxes[k].xyz);
        float shape = uQuakeAxes[k].w < 2.5 ? 0.5 * (3.0 * x * x - 1.0) : 0.5 * x * (5.0 * x * x - 3.0);
        f *= 1.0 + 2.0 * swing * shape;
    }
    return max(f, 0.3);
}
#endif

// A direction in Kerr-Schild Cartesian space, carried into the body's rotating frame
// (where its surface pattern is painted) as it was `delay` ago.
vec3 toRotatingFrame(vec3 d, float delay, float r) {
    float turn = uSurfacePhase - uSurfaceOmega * delay + uSurfaceShift + atan(spinA, r);
    float c = cos(turn), s = sin(turn);
    return vec3(c * d.x + s * d.y, -s * d.x + c * d.y, d.z);
}

// Where the chord X0 -> X1 first meets the surface r = R, as a fraction along it (2 if
// it does not). In Kerr-Schild Cartesian space r = R is the spheroid
// (x^2 + y^2) / (R^2 + a^2) + z^2 / R^2 = 1: scaled to a sphere, an exact line test, so
// a ray grazing the limb between two steps still finds the surface.
float surfaceCrossing(vec3 X0, vec3 X1, float r1) {
    float R = uSurfaceR;
    float squash = R * inversesqrt(R * R + spinA2);
    vec3 k = vec3(squash, squash, 1.0);
    vec3 a = X0 * k, b = (X1 - X0) * k;
    float bb = dot(b, b), ab = dot(a, b), c = dot(a, a) - R * R;
    float disc = ab * ab - bb * c;
    if (bb > 0.0 && disc >= 0.0) {
        float t = (-ab - sqrt(disc)) / bb;
        if (t >= 0.0 && t <= 1.0) return t;
    }
    return r1 <= R ? 1.0 : 2.0;
}

// --- The air of an Earth-like planet ------------------------------------------------
//
// Single scattering of sunlight: Rayleigh scattering by molecules (scale height 8 km,
// 1/lambda^4: the blue sky and red sunsets), Mie scattering by aerosols (1.2 km, mostly
// forward), and absorption by ozone (the Chappuis band, which keeps the twilight sky
// blue). Earth's coefficients at sea level, in 1/m (Bruneton and Neyret 2008; Hillaire
// 2020). Ozone really lies in a layer near 25 km; here it follows the Rayleigh profile
// with the same column (optical depth 0.010, 0.028, 0.001). Drawn at true scale: the
// air is 1.3% of the Earth's radius thick (100 km) and its density halves every 5.5 km.
const vec3 kRayleighSI = vec3(5.802e-6, 13.558e-6, 33.1e-6);
const vec3 kOzoneSI = vec3(1.22e-6, 3.53e-6, 0.16e-6);
const float kMieScatterSI = 3.996e-6, kMieExtinctSI = 4.40e-6;
const float kRayleighHeightSI = 8000.0, kMieHeightSI = 1200.0;
const float kMieG = 0.8;
const float kAirTopSI = 100000.0;

// The Chapman function times exp(-h): optical depth from altitude h along zenith angle
// chi through an exponential atmosphere on a sphere, in units of the vertical column
// from the ground (Schuler 2012, "An approximation to the Chapman grazing-incidence
// function for atmospheric scattering"). X = R / H50 and h = altitude / H50, where H50
// = H ln 2 is the height over which the density halves.
float chapmanH(float X, float h, float coschi) {
    float c = sqrt(X + h);
    if (coschi >= 0.0) return c / (c * coschi + 1.0) * exp2(-h);
    float x0 = sqrt(1.0 - coschi * coschi) * (X + h);
    float c0 = sqrt(x0);
    return 2.0 * c0 * exp2(X - x0) - c / (1.0 - c * coschi) * exp2(-h);
}

// Transmittance of the air from a point at radius r (render units) toward the sun:
// zero in the planet's shadow.
vec3 airToSun(float r, float cosChi, float metersPerUnit) {
    float R = uSurfaceR;
    // Below the geometric horizon the sunlight passes through the planet.
    float horizon = -sqrt(max(1.0 - (R * R) / (r * r), 0.0));
    float lit = smoothstep(horizon - 0.002, horizon + 0.002, cosChi);
    if (lit <= 0.0) return vec3(0.0);
    float hR = kRayleighHeightSI / metersPerUnit, hM = kMieHeightSI / metersPerUnit;
    float alt = max(r - R, 0.0);
    float X50R = R / (hR * 0.6931472), X50M = R / (hM * 0.6931472);
    float chR = chapmanH(X50R, alt / (hR * 0.6931472), max(cosChi, horizon));
    float chM = chapmanH(X50M, alt / (hM * 0.6931472), max(cosChi, horizon));
    vec3 tau = (kRayleighSI + kOzoneSI) * kRayleighHeightSI * chR + vec3(kMieExtinctSI * kMieHeightSI * chM);
    return exp(-tau) * lit;
}

// Sunlight scattered toward the camera by the air along the line o + t d, between t0 and
// t1 (around a planet, light bends by ~1e-9 rad: a straight line). `trans` returns the
// transmittance of the segment, which dims whatever lies behind it. sun: the sun's
// color as it lights a white surface (radiance of a white Lambertian surface it faces).
vec3 airGlow(vec3 o, vec3 d, float t0, float t1, vec3 sun, float jitter, out vec3 trans) {
    trans = vec3(1.0);
    if (t1 <= t0) return vec3(0.0);
    float R = uSurfaceR;
    float metersPerUnit = uSurfaceKm * 1000.0 / R;
    vec3 betaR = kRayleighSI * metersPerUnit, betaO = kOzoneSI * metersPerUnit;
    float betaMs = kMieScatterSI * metersPerUnit, betaMe = kMieExtinctSI * metersPerUnit;
    float hR = kRayleighHeightSI / metersPerUnit, hM = kMieHeightSI / metersPerUnit;
    vec3 L = uSunDirection;
    float mu = dot(d, L);  // cosine of the scattering angle
    float phaseR = 3.0 / (16.0 * PI) * (1.0 + mu * mu);
    float g2 = kMieG * kMieG;
    float phaseM = (1.0 - g2) / (4.0 * PI * pow(1.0 + g2 - 2.0 * kMieG * mu, 1.5));
    // The densest air is at the line's lowest point: split there and crowd the samples
    // toward it on both sides.
    float tc = clamp(-dot(o, d), t0, t1);
    vec3 inscatter = vec3(0.0);
    const int N = 8;
    for (int part = 0; part < 2; ++part) {
        float span = part == 0 ? tc - t0 : t1 - tc;
        if (span <= 0.0) continue;
        for (int i = 0; i < N; ++i) {
            float v = (float(i) + jitter) / float(N);
            float t, dt;
            if (part == 0) {
                t = t0 + span * (1.0 - (1.0 - v) * (1.0 - v));
                dt = span * 2.0 * (1.0 - v) / float(N);
            } else {
                t = tc + span * v * v;
                dt = span * 2.0 * v / float(N);
            }
            vec3 x = o + t * d;
            float r = length(x);
            float alt = r - R;
            float rhoR = exp(-alt / hR), rhoM = exp(-alt / hM);
            vec3 extinction = (betaR + betaO) * rhoR + vec3(betaMe * rhoM);
            vec3 stepTrans = exp(-extinction * dt);
            vec3 sunT = airToSun(r, dot(x, L) / r, metersPerUnit);
            vec3 scatter = (betaR * rhoR * phaseR + vec3(betaMs * rhoM * phaseM)) * sunT;
            // Exact for constant coefficients over the step: (1 - e^-sigma dt) / sigma.
            inscatter += trans * scatter * (1.0 - stepTrans) / max(extinction, vec3(1e-12));
            trans *= stepTrans;
        }
    }
    return inscatter * sun * PI;
}

// Where the line o + t d crosses a sphere of radius rad: (enter, leave); enter > leave
// if it misses.
vec2 sphereSpan(vec3 o, vec3 d, float rad) {
    float b = dot(o, d), c = dot(o, o) - rad * rad;
    float disc = b * b - c;
    if (disc < 0.0) return vec2(1.0, -1.0);
    float s = sqrt(disc);
    return vec2(-b - s, -b + s);
}

// --- The solar corona ---------------------------------------------------------------
//
// The white-light (K) corona is photospheric light Thomson-scattered by free electrons.
// Electron density from Baumbach and Allen, n_e = 1e8 (2.99 rho^-16 + 1.55 rho^-6 +
// 0.036 rho^-1.5) per cm^3 (rho = r / R), times the dilution of the photosphere's light
// W = (1 - sqrt(1 - 1/rho^2)) / 2; integrated along the line of sight this reproduces
// Baumbach's brightness profile within 5% from 1.02 to 6 R (src/surface.cpp, tested).
// The magnetic field shapes it: helmet streamers over the equator and the active
// regions, bulbous low down and narrowing outward into thin rays; plumes over the
// poles. uCoronaScale converts to the disk center's brightness (1e-6 of it at 1.1 R for
// the real Sun, times the chosen exaggeration).
float electronDensity(float rho) {
    float l = log(rho);
    return 2.99 * exp(-16.0 * l) + 1.55 * exp(-6.0 * l) + 0.036 * exp(-1.5 * l);
}

float coronaStructure(vec3 dir, float rho, float fp) {
    float lat = asin(clamp(dir.z, -1.0, 1.0));
    float width = 0.32 / rho + 0.04;  // streamers narrow outward
    float wave = 0.15 * gradientNoise(vec3(dir.xy * 1.5, 3.0));
    float eq = (lat - wave) / width;
    float streamers = exp(-eq * eq) * (0.55 + 0.45 * gradientNoise(dir * 2.5 + 11.0));
    float overRegions = 0.0;
    for (int k = 0; k < uRegionCount; ++k) {
        vec3 rd = uRegions[k].xyz;
        float a = sqrt(max(2.0 - 2.0 * dot(dir, rd), 0.0)) / (width * 0.9 + uRegions[k].w);
        overRegions += exp(-a * a);
    }
    float plumes = smoothstep(0.82, 0.95, abs(dir.z)) * (0.5 + 0.5 * gradientNoise(dir * 45.0));
    // Fine radial striation: brightness that depends on direction alone runs out
    // along the radius.
    float rays = 1.0 + 0.35 * fbmFiltered(dir * 18.0, 3, fp * 18.0);
    return (0.22 + 1.7 * streamers + 1.1 * overRegions + 0.5 * plumes) * rays;
}

// The corona's light along the line o + t d between t0 and t1 (a straight line: around
// an ordinary star light bends by microradians). In units of the disk center's
// radiance. The density falls so steeply that nearly all the light comes from near the
// line's closest approach to the star, so the streamers' pattern is read there, once.
float coronaGlow(vec3 o, vec3 d, float t0, float t1, float jitter, float delay) {
    if (t1 <= t0) return 0.0;
    float R = uSurfaceR;
    float tc = clamp(-dot(o, d), t0, t1);
    float sum = 0.0;
    const int N = 8;
    for (int part = 0; part < 2; ++part) {
        float span = part == 0 ? tc - t0 : t1 - tc;
        if (span <= 0.0) continue;
        for (int i = 0; i < N; ++i) {
            float v = (float(i) + jitter) / float(N);
            float t, dt;
            if (part == 0) {
                t = t0 + span * (1.0 - (1.0 - v) * (1.0 - v));
                dt = span * 2.0 * (1.0 - v) / float(N);
            } else {
                t = tc + span * v * v;
                dt = span * 2.0 * v / float(N);
            }
            float rho = max(length(o + t * d) / R, 1.0001);
            float W = 0.5 * (1.0 - sqrt(max(1.0 - 1.0 / (rho * rho), 0.0)));
            sum += electronDensity(rho) * W * dt / R;
        }
    }
    vec3 xc = o + tc * d;
    float rc = length(xc);
    float rhoC = max(rc / R, 1.0001);
    // A pixel's size there, as an angle seen from the star's center.
    float fpDir = tc * uPixelAngle / rc;
    return sum * coronaStructure(toRotatingFrame(xc / rc, delay, rc), rhoC, fpDir) * uCoronaScale;
}

// The chromosphere: a thin layer (~2,000 km on the Sun) above the photosphere that
// glows in the red H-alpha line and blue H-beta: the pink rim seen in eclipses, made
// ragged by spicules. Seen against the disk it is lost in the photosphere's glare;
// beyond the limb, a grazing line through a thin exponential layer collects
// exp(-h/H) sqrt(2 pi R H) (h: the line's lowest altitude).
vec3 chromosphere(vec3 o, vec3 d, float delay) {
    float R = uSurfaceR;
    float tc = -dot(o, d);
    if (tc <= 0.0) return vec3(0.0);
    vec3 x = o + tc * d;
    float r = length(x);
    float h = r - R;
    if (h < 0.0) return vec3(0.0);
    vec3 dir = toRotatingFrame(x / r, delay, r);
    float H = uChromosphereHeight * R * (1.0 + 0.6 * gradientNoise(dir * 260.0) + 0.3 * gradientNoise(dir * 60.0));
    // In film colors, the spicule forest of ultraviolet images: taller and more ragged.
    if (uFilmColor == 1) H *= 2.5 * (0.6 + 0.8 * smoothstep(-0.2, 0.5, gradientNoise(dir * 420.0)));
    H = max(H, 1e-6 * R);
    float column = exp(-h / H) * sqrt(2.0 * PI * R * H) / R;
    return vec3(1.0, 0.16, 0.3) * column * uChromosphereScale;
}

#if SURFACE_LOOK == 1
#include "sdo.glsl"

// The baked 30.4 nm map (sdo.glsl), equirectangular in the rotating frame.
layout(binding = 11) uniform sampler2D uSdoMap;
float sdoSample(vec3 n, float fp) {
    float u = 0.5 + atan(n.y, n.x) / (2.0 * PI);
    float v = 0.5 + asin(clamp(n.z, -1.0, 1.0)) / PI;
    float lod = log2(max(fp * float(textureSize(uSdoMap, 0).x) / (2.0 * PI), 1.0));
    return textureLod(uSdoMap, vec2(u, v), lod).r;
}

// --- A red giant, as ESO's images and simulations of Betelgeuse show it -----------------------
//
// A red supergiant's surface holds only a handful of convection cells, each a good part
// of the star across (Schwarzschild 1975; the radiation-hydrodynamics models of Freytag
// et al.): broad, soft, bright plumes of rising gas separated by darker, cooler lanes, on
// top of a few even larger hot regions, the whole slowly churning (at the cells'
// lifetime, months). The limb is darkened deeply and blurred by an extended atmosphere
// (giantEnvelope). Intensity 1 is a bright plume.
float giantIntensity(vec3 n, float mu, float fp, float phase) {
    // Interferometry and the models show fewer, larger cells than the pressure scale
    // height alone suggests: a handful across the disk.
    float freq = clamp(0.3 * uGranuleFreq, 1.5, 4.0);
    vec3 p = n * freq;
    vec3 warp = vec3(gradientNoise(p * 0.6 + 3.0), gradientNoise(p * 0.6 + 8.0), gradientNoise(p * 0.6 + 1.0));
    vec3 w = worley(p + 0.9 * warp, phase * 0.15);
    float plumes = 1.0 - smoothstep(0.0, 0.75, w.x);             // bright, domed cell centers
    float lanes = 1.0 - smoothstep(0.0, 0.4, w.y - w.x);         // the cool gaps between them
    float big = smoothstep(-0.35, 0.55, fbmFiltered(n * 1.6 + vec3(5.0, 1.0, phase * 0.03), 3, fp * 1.6));
    float fine = fbmFiltered(p * 5.0 + 2.0, 5, fp * freq * 5.0);
    float I = 0.42 + 0.45 * plumes * (0.55 + 0.45 * big) + 0.25 * big - 0.18 * lanes * lanes + 0.14 * fine;
    return I * (0.18 + 0.82 * pow(mu, 0.55));
}

// Light from a red giant's extended atmosphere (its molecular layers and the dust it is
// shedding) along the line o + t d (t0 to t1): a glow falling off steeply over the first
// ~0.2 radii and then slowly, broken into plumes of gas thrown off at a few spots (ESO's
// images of Betelgeuse show such a plume, six times the star's radius). In units of a
// bright plume's brightness.
float giantEnvelope(vec3 o, vec3 d, float t0, float t1, float jitter) {
    if (t1 <= t0) return 0.0;
    float R = uSurfaceR;
    float turn = uSurfacePhase + uSurfaceShift;
    float ct = cos(turn), st = sin(turn);
    float sum = 0.0;
    const int N = 6;
    float tc = clamp(-dot(o, d), t0, t1);
    for (int part = 0; part < 2; ++part) {
        float span = part == 0 ? tc - t0 : t1 - tc;
        if (span <= 0.0) continue;
        for (int i = 0; i < N; ++i) {
            float v = (float(i) + jitter) / float(N);
            float t, dt;
            if (part == 0) {
                t = t0 + span * (1.0 - (1.0 - v) * (1.0 - v));
                dt = span * 2.0 * (1.0 - v) / float(N);
            } else {
                t = tc + span * v * v;
                dt = span * 2.0 * v / float(N);
            }
            vec3 x = o + t * d;
            float r = length(x);
            float h = max(r / R - 1.0, 0.0);
            vec3 dir = x / r;
            vec3 rot = vec3(ct * dir.x + st * dir.y, -st * dir.x + ct * dir.y, dir.z);
            // Clumpy gas, streaming out (three octaves: this is read many times per ray).
            vec3 q = rot * 3.0 + vec3(0.0, 0.0, -h * 2.0);
            float noise = (0.5 * valueNoise(q) + 0.25 * valueNoise(q * 2.07 + 13.7) + 0.125 * valueNoise(q * 4.3 + 7.3)) / 0.875;
            float clumps = 0.4 + 0.9 * smoothstep(0.2, 0.75, noise);
            float plume = 0.0;
            for (int k = 0; k < 2; ++k) {
                vec3 pd = normalize(vec3(cos(1.3 + 3.4 * float(k)), sin(1.3 + 3.4 * float(k)), 0.55 - 0.9 * float(k)));
                float a = 1.0 - dot(rot, pd);
                plume += exp(-a / (0.015 + 0.02 * h)) * exp(-h / 0.7);
            }
            float density = exp(-h / 0.07) + 0.15 * exp(-h / 0.35) * clumps + 0.18 * plume * clumps;
            sum += density * dt / R;
        }
    }
    return sum;
}

// --- Prominences and coronal loops -------------------------------------------------------
//
// Over an active region the field rises out of one polarity of the spot group and arches
// back into the other: an arcade of loops side by side, along the group (tilted toward
// the equator by half its latitude, Joy's law), holding gas at ~10,000 K that glows in
// H-alpha and the ultraviolet (the loops and prominences of eclipses and of SDO's 30.4 nm
// images). Each loop is an arc of a circle standing on the surface, ~5,000 km thick and
// up to a tenth of the radius high. Quiescent prominences, long sheets held up in dips of
// the field, hang over the high-latitude filament channels (the polar crown).
// Each lies in a plane, so a straight line of sight crosses it once: the glow is read
// there, from the tube's (or sheet's) profile, exactly whatever the pixel size. Seen
// against the disk they are lost in its glare (in H-alpha they would show as dark
// filaments).
const int kLoopsPerArcade = 4;
const float kLoopWidth = 0.0035;  // half-width of a loop, in surface radii (~2,500 km on the Sun)

// The glow of a tube of Gaussian half-width w around the circle (center c, radius a) in
// the plane of b and n (normal m), collected by the line o + t d near the point at t.
// Locally the tube is straight: its profile integrates across to sqrt(pi) w / sin(gamma)
// exp(-D^2 / w^2), gamma the angle between line and tube, D their closest distance.
// Farther than sqrt(2 a w) along the line the tube has curved away by w: there the
// distance at t is used instead. pc: the circle's point nearest; 0 if it is underground.
float tubeGlow(vec3 o, vec3 d, float t, float tEnd, vec3 c, vec3 b, vec3 n, float a, float w, out vec3 pc,
               out vec2 across) {
    across = vec2(0.0);
    vec3 x = o + t * d;
    vec3 q = x - c;
    float u = dot(q, b), v = dot(q, n);
    float rho = max(sqrt(u * u + v * v), 1e-6);
    pc = c + (a / rho) * (u * b + v * n);
    if (dot(pc, pc) < uSurfaceR * uSurfaceR) return 0.0;
    vec3 T = (-v * b + u * n) / rho;  // along the tube
    vec3 delta = x - pc;
    float dT = dot(d, T);
    float sin2 = max(1.0 - dT * dT, 0.02);
    float s = (dot(delta, T) * dT - dot(delta, d)) / sin2;  // to the closest approach
    if (t + s < 0.0 || t + s > tEnd) return 0.0;
    vec3 gap = delta + s * d;
    gap -= dot(gap, T) * T;
    across = vec2(dot(gap, pc - c) / a, dot(gap, cross(b, n))) / w;
    float D2 = mix(dot(gap, gap), dot(delta, delta), smoothstep(0.5, 1.5, abs(s) / sqrt(2.0 * a * w)));
    return 1.7725 * w * inversesqrt(sin2) * exp(-D2 / (w * w));
}

float prominences(vec3 o, vec3 d, float tEnd, bool overDisk) {
    float R = uSurfaceR;
    float turn = uSurfacePhase + uSurfaceShift + atan(spinA, R);
    float ct = cos(turn), st = sin(turn);
    float glow = 0.0;
    float time = uGranulePhase * 0.02;  // the gas drains and refills over many granule lifetimes
    for (int k = 0; k < uRegionCount + kQuiescentProminences; ++k) {
        vec3 nRot;
        float size;
        bool quiescent = k >= uRegionCount;
        bool sheet = false;
        if (!quiescent) {
            nRot = uRegions[k].xyz;
            size = 1.3 * uRegions[k].w;
        } else {
            vec4 q = kQuiescent[k - uRegionCount];
            nRot = vec3(cos(q.x) * cos(q.y), cos(q.x) * sin(q.y), sin(q.x));
            size = q.z;
            sheet = q.w > 0.5;
        }
        // Into the inertial frame (the star has turned through `turn`).
        vec3 n = vec3(ct * nRot.x - st * nRot.y, st * nRot.x + ct * nRot.y, nRot.z);
        // Quick reject: the line of sight passes nowhere near this region's loops.
        vec3 toRegion = n * R - o;
        float along = dot(toRegion, d);
        float miss = length(toRegion - along * d);
        if (miss > 1.6 * size * R + 0.02 * R) continue;
        vec3 east = cross(vec3(0.0, 0.0, 1.0), n);
        east = dot(east, east) > 1e-8 ? normalize(east) : vec3(1.0, 0.0, 0.0);
        vec3 north = cross(n, east);
        float lat = asin(clamp(n.z, -1.0, 1.0));
        float tilt = sheet ? 0.25 * lat : 0.5 * lat;
        vec3 b = east * cos(tilt) - north * sin(tilt);  // along the arcade's loops
        vec3 m = cross(n, b);                           // the loops' planes are normal to this
        float dm = dot(d, m);
        int loops = quiescent ? 1 : kLoopsPerArcade;
        for (int i = 0; i < loops; ++i) {
            vec4 h = hash43(vec3(float(k), float(i), 5.1));
            float a = size * R * (quiescent ? 0.9 : 0.4 + 0.35 * h.x);
            // The circle's center sinks half its radius under the surface: an arch half as
            // high as its footpoints are apart, as loops and prominences stand.
            vec3 c = n * (R - 0.5 * a) + m * ((float(i) - 0.5 * float(loops - 1)) * 0.22 * size * R);
            vec3 q0 = o - c;
            if (sheet) {
                // A hedgerow: a sheet filling the arch, in vertical threads, ragged on top,
                // read where the line pierces its plane.
                if (abs(dm) < 1e-5) continue;
                float t = -dot(q0, m) / dm;
                if (t < 0.0 || t > tEnd) continue;
                vec3 x = o + t * d;
                float r = length(x);
                if (r < R) continue;
                vec3 q = x - c;
                float rho = length(q);
                float alpha = atan(dot(q, n), dot(q, b));
                float w = max(kLoopWidth * R * 2.0, 0.7 * uPixelAngle * t);
                float top = a * (0.85 + 0.15 * gradientNoise(vec3(alpha * 5.0, h.z * 10.0, time)));
                float inside = smoothstep(top + w, top - w - 0.15 * a, rho);
                float threads = 0.45 + 0.55 * smoothstep(-0.3, 0.6, gradientNoise(vec3(dot(q, b) / (0.012 * R), h.w * 7.0, time)));
                float height = (r - R) / R;
                glow += inside * threads * 0.45 / max(abs(dm), 0.25) * exp(-height / 0.3) * (0.6 + 0.8 * h.w);
                continue;
            }
            // A loop: a thin tube along the arc. A line crossing its plane steeply meets it
            // where it pierces the plane; one running along the plane, where its projection
            // crosses the circle (or passes closest to it).
            float w0 = kLoopWidth * R * (quiescent ? 6.0 : 0.7 + 0.6 * h.y);
            float w = max(w0, 0.7 * uPixelAngle * (dot(c - o, d)));
            float steep = smoothstep(0.2, 0.4, abs(dm));
            float candidates[3] = float[3](-1.0, -1.0, -1.0);
            float weights[3] = float[3](0.0, 0.0, 0.0);
            if (steep > 0.0) {
                candidates[0] = -dot(q0, m) / dm;
                weights[0] = steep;
            }
            if (steep < 1.0) {
                float du = dot(d, b), dv = dot(d, n);
                float u0 = dot(q0, b), v0 = dot(q0, n);
                float P = du * du + dv * dv, Q = u0 * du + v0 * dv, S = u0 * u0 + v0 * v0 - a * a;
                float disc = Q * Q - P * S;
                if (P > 1e-8) {
                    if (disc >= 0.0) {
                        float sq = sqrt(disc);
                        candidates[1] = (-Q - sq) / P;
                        candidates[2] = (-Q + sq) / P;
                        weights[1] = weights[2] = 1.0 - steep;
                    } else {
                        candidates[1] = -Q / P;
                        weights[1] = 1.0 - steep;
                    }
                }
            }
            for (int j = 0; j < 3; ++j) {
                if (weights[j] <= 0.0) continue;
                vec3 pc;
                vec2 across;
                float tube = tubeGlow(o, d, candidates[j], tEnd, c, b, n, a, w, pc, across);
                if (tube <= 0.0) continue;
                // Threads and knots of gas strung along the arc; brightest low down, where
                // the gas is densest.
                vec3 q = pc - c;
                float alpha = atan(dot(q, n), dot(q, b));
                float knots = 0.35 + 0.65 * smoothstep(-0.4, 0.5, gradientNoise(vec3(alpha * 4.0, h.z * 10.0, time)));
                // A big prominence is a bundle of threads along the field.
                if (quiescent) knots *= 0.45 + 0.9 * smoothstep(-0.35, 0.45, gradientNoise(vec3(alpha * 2.5, across * 1.6 + h.xy * 9.0)));
                float height = (length(pc) - R) / R;
                // In units of the peak glow of a tube seen across.
                glow += weights[j] * tube / (1.7725 * w) * (w0 / w) * 2.0 * knots *
                        exp(-height / (quiescent ? 0.3 : 0.08)) * (0.6 + 0.8 * h.w);
            }
        }
    }
    if (overDisk) glow *= 0.12;
    return 0.6 * glow;
}

#endif

#if SURFACE_LOOK == 5 || SURFACE_LOOK == 6
// --- Real maps -------------------------------------------------------------------------------
//
// When found (assets/planets), the planets are drawn from spacecraft maps: for the Earth,
// NASA's Blue Marble Next Generation (December 2004, with the sea floor), its cloud cover
// and the Black Marble's city lights (2016); for Jupiter, Cassini's global map (PIA07782,
// December 2000). Equirectangular, east to the right, north up; read at the mipmap level
// of the pixel's size on the surface. Everything else (the light, the air, the glint, the
// shadows, the night side) is still computed.
layout(binding = 8) uniform sampler2D uMapDay;      // sRGB: linear when read
layout(binding = 9) uniform sampler2D uMapClouds;   // cover, 0..1
layout(binding = 10) uniform sampler2D uMapNight;   // city lights
uniform int uPlanetMaps;                            // 1: the maps are loaded

vec4 mapSample(sampler2D map, vec3 n, float fp) {
    float u = 0.5 + atan(n.y, n.x) / (2.0 * PI);
    float v = 0.5 + asin(clamp(n.z, -1.0, 1.0)) / PI;
    float lod = log2(max(fp * float(textureSize(map, 0).x) / (2.0 * PI), 1.0));
    return textureLod(map, vec2(u, v), lod);
}
float mapClouds(vec3 n, float fp) {
    return 0.97 * pow(mapSample(uMapClouds, n, fp).r, 1.4);
}
#endif

// --- The surface ----------------------------------------------------------------------

// Light leaving the surface toward the camera, `delay` ago, at (R, theta, phi). hitX is
// the point in Kerr-Schild Cartesian space, rayDir the backward ray's direction there,
// and dist how far the light travelled to the camera (it sets the size of a pixel on
// the surface).
vec3 surfaceRadiance(float th, float ph, float delay, vec3 rayDir, vec3 hitX, float dist) {
    float R = uSurfaceR;
    // The surface element rotates rigidly: u = u^t (1, 0, 0, Omega), and the light
    // reaching us was emitted at g = 1 / (u^t (E - Omega L)) (as for the disk).
    float gtt, gtp, gpp;
    blMetric(R, th, gtt, gtp, gpp);
    float ut = inversesqrt(max(-(gtt + 2.0 * uSurfaceOmega * gtp + uSurfaceOmega * uSurfaceOmega * gpp), 1e-12));
    float g = 1.0 / max(ut * (photonE - uSurfaceOmega * photonL), 1e-4);
    // Emission angle from the normal, in the surface's own frame:
    // cos = (p . e_r) / (-p . u) = g |v_r| / sqrt(Sigma Delta).
    float c = cos(th), s = sin(th);
    float sigma = R * R + spinA2 * c * c;
    float mu = clamp(g * sqrt(max(radialPotential(R), 0.0)) / sqrt(max(sigma * deltaAt(R), 1e-12)), 0.0, 1.0);
    // Where on the body, in its rotating frame: longitude from phi_BL, less the angle it
    // had turned through when the light left.
    float turned = uSurfacePhase - uSurfaceOmega * delay;
    float lon = ph - uSurfaceShift - turned;
    vec3 n = vec3(s * cos(lon), s * sin(lon), c);
    float flash = uFlash.x - delay >= 0.0 ? uFlash.y * exp(-(uFlash.x - delay) / uFlash.z) : 0.0;
    // A pixel's size on the surface, in radians, stretched toward the limb.
    float fp = dist * uPixelAngle / R / max(mu, 0.12);

#if SURFACE_LOOK == 5 || SURFACE_LOOK == 6
    {
        // A planet shines by reflected light, from a distant star (and faintly a disk).
        vec3 normal = normalize(hitX);
        vec3 V = -rayDir;
        vec3 L = uSunDirection;
        float mu0 = dot(normal, L);
        float muV = max(dot(normal, V), 0.0);
        vec3 sun = shiftedBlackbody(5772.0, g, 5772.0) * uSurfaceBrightness * uPlanetLight;
        vec3 diskLight = uDiskEnabled == 1 ? vec3(0.15 * uDiskBrightness * abs(normal.z)) * sun : vec3(0.0);
#if SURFACE_LOOK == 6
        vec3 albedo = jupiterAlbedo(n, fp, uPlanetSeed);
        if (uPlanetMaps == 1) {
            // Cassini's colors as reflectance (Jupiter's albedo is ~0.5); past ~78 degrees,
            // where the map has no data, the drawn polar regions take over.
            vec3 map = mapSample(uMapDay, n, fp).rgb * 0.85;
            albedo = mix(map, albedo, smoothstep(1.32, 1.42, abs(asin(clamp(n.z, -1.0, 1.0)))));
        }
        vec3 lit = sun * minnaert(mu0, muV, 0.9) * smoothstep(-0.03, 0.05, mu0) + diskLight + 0.002 * sun;
        return albedo * lit * (1.0 + flash);
#else
        vec3 Lr = toRotatingFrame(L, delay, R);
        Ground gr;
        if (uPlanetMaps == 1) {
            // The Blue Marble as surface reflectance (the image is brightened for display:
            // deserts ~0.35, as measured), its seas found by their color for the glint;
            // NASA's cloud cover; the Black Marble's lights above its moonlit ground.
            vec3 day = mapSample(uMapDay, n, fp).rgb;
            gr.albedo = day * 0.6;
            gr.water = smoothstep(0.004, 0.02, day.b - 1.2 * max(day.r, day.g)) *
                       (1.0 - smoothstep(0.06, 0.12, dot(day, vec3(0.2126, 0.7152, 0.0722))));
            gr.clouds = mapClouds(n, fp);
            gr.lights = 1.6 * max(dot(mapSample(uMapNight, n, fp).rgb, vec3(0.3, 0.5, 0.2)) - 0.16, 0.0);
        } else {
            gr = earthGround(n, fp, uPlanetSeed, Lr);
        }
        float metersPerUnit = uSurfaceKm * 1000.0 / R;
        vec3 sunT = uAtmosphere == 1 ? airToSun(R, mu0, metersPerUnit) : vec3(step(0.0, mu0));
        // Light from the blue sky fills in the shadows; it fades through twilight.
        vec3 skyFill = uAtmosphere == 1 ? vec3(0.025, 0.05, 0.11) * smoothstep(-0.12, 0.35, mu0) : vec3(0.0);
        // Clouds cast shadows: is there cloud ~6 km up, toward the sun?
        float shade = 1.0;
        if (mu0 > 0.0 && gr.clouds < 0.98) {
            vec3 toward = normalize(n + Lr * (6.0e3 / (uSurfaceKm * 1000.0)) / max(mu0, 0.05));
            shade = 1.0 - 0.8 * (uPlanetMaps == 1 ? mapClouds(toward, fp) : earthClouds(toward, fp, uPlanetSeed));
        }
        vec3 direct = sun * sunT * max(mu0, 0.0);
        vec3 groundLight = gr.albedo * (direct * shade + sun * skyFill + diskLight);
        groundLight += gr.water * (1.0 - gr.clouds) * shade * sunGlint(normal, L, V, 0.17) * sun * sunT;
        vec3 cloudLight = vec3(0.86) * (direct + sun * skyFill + diskLight);
        vec3 color = mix(groundLight, cloudLight, gr.clouds);
        // Cities at night, sodium-orange (brightened, as in long exposures: real ones are
        // ~1e-5 of the day side).
        float night = (1.0 - smoothstep(-0.12, 0.08, mu0));
        color += night * gr.lights * (1.0 - 0.75 * gr.clouds) * vec3(1.0, 0.6, 0.27) * 0.6 * uSurfaceBrightness;
        return color * (1.0 + flash);
#endif
    }
#endif

    // Thermal light. Temperature factors are applied to the physical temperature; above
    // 30,000 K the display temperature goes as T^0.3 (central::displayTemperature).
    float tFactor = 1.0;
    float limbBoost = 1.0;  // faculae
    float filmIntensity = 1.0;  // brightness in film colors (sdoIntensity)
    if (uSurfaceDetail == 1) {
#if SURFACE_LOOK == 1
        {
            float phase = uGranulePhase - uGranuleRate * delay;
            tFactor *= 1.0 + uGranuleContrast * granulation(n, uGranuleFreq, phase, fp, uGranuleLane);
            float network = supergranuleNetwork(n, uGranuleFreq / 25.0, fp);
            float penumbra;
            tFactor *= spotTemperature(n, uSpotDrop, fp, penumbra);
            float plage = faculae(n, network);
            limbBoost = 1.0 + (1.0 - mu) * (0.06 * network + 0.12 * plage) * (1.0 - penumbra);
            if (uFilmColor == 1) {
                if (uGranuleLane > 0.3) {
                    filmIntensity = giantIntensity(n, mu, fp, phase) * pow(tFactor, 2.0);
                } else {
                    // The limb brightens: there the line of sight runs through more of the thin layer.
                    filmIntensity = sdoSample(n, fp) * (1.0 + 0.5 * pow(1.0 - mu, 3.0)) * mix(1.0, pow(tFactor, 4.0), 0.5);
                }
            }
        }
#elif SURFACE_LOOK == 2
        tFactor *= brownDwarfWeather(n, fp, uPlanetSeed);
#elif SURFACE_LOOK == 4
        {
            // Heat flows along the magnetic field far more easily than across it, so a
            // dipole's poles are hotter than its equator (Greenstein and Hartke 1983):
            // T^4 = Tp^4 cos^2(Theta_B) + Teq^4 sin^2(Theta_B), Theta_B the angle between
            // field and normal (cos^2 = 4 cos^2 / (1 + 3 cos^2) of the magnetic
            // colatitude), Teq = 0.6 Tp; normalized to the mean surface temperature.
            float cm = dot(n, uMagneticAxis);
            float cos2B = 4.0 * cm * cm / (1.0 + 3.0 * cm * cm);
            tFactor *= pow((0.13 + 0.87 * cos2B) / 0.5885, 0.25);
            if (uQuakeCount > 0) tFactor *= quakeShimmer(n, delay);
            if (uMagnetar == 1) {
                // A magnetar's twisted field drives currents that heat one patch of the
                // crust: X-ray spectra show a hot spot ~twice the surface temperature
                // covering several percent of it (as on XTE J1810-197), off the dipole
                // pole. The rest is smooth: a neutron star's crust has no relief to speak of.
                vec3 spotAxis = normalize(uMagneticAxis + vec3(0.35, -0.2, 0.0));
                float spot = smoothstep(0.90, 0.95, dot(n, spotAxis));
                tFactor *= 1.0 + 0.8 * spot;
            }
        }
#endif
    }
    float tExp = uSurfaceTemp > 30000.0 ? 0.3 : 1.0;
    float T = uSurfaceTemp * pow(tFactor * limbTemperature(mu), tExp);
    float Tref = uSurfaceTemp * pow(kCenterTemperature, tExp);
    vec3 color = shiftedBlackbody(T, g, Tref) * (kPhotosphereLevel * uSurfaceBrightness * limbBoost);
#if SURFACE_LOOK == 1
    if (uFilmColor == 1) {
        // Cooler stars redder in the same table; stars hotter than ~7,000 K keep the
        // blackbody palette.
        float shift = clamp(pow(uSurfaceTemp / 5772.0, 0.6), 0.7, 1.0);
        if (uGranuleLane > 0.3) shift = 1.25;  // giantIntensity is scaled to its own range
        vec3 sdo = aia304(filmIntensity * 0.43 * shift * (0.5 + 0.5 * g)) * (kSdoLevel * uSurfaceBrightness);
        color = mix(sdo, filmBlackbody(T, g, Tref) * (uSurfaceBrightness * limbBoost),
                    smoothstep(6500.0, 8500.0, uSurfaceTemp));
    }
#elif SURFACE_LOOK == 2
    if (uFilmColor == 1) {
        // As NASA's and films' brown dwarfs: a dark maroon ball wrapped in glowing orange
        // and pink cloud streaks, thin and wind-sheared along the latitude bands, through
        // whose gaps the hotter depths shine (the bands of the weather above, sharpened).
        float lat = asin(clamp(n.z, -1.0, 1.0));
        float lon = atan(n.y, n.x);
        float shear = fbmFiltered(n * vec3(2.5, 2.5, 7.0) + uPlanetSeed, 5, fp * 7.0);
        float streaks = ridgedFiltered(vec3(cos(lon) * 1.5, sin(lon) * 1.5, lat * 22.0 + 3.0 * shear) + uPlanetSeed, 5,
                                       fp * 22.0);
        float wisps = fbmFiltered(vec3(cos(lon) * 4.0, sin(lon) * 4.0, lat * 40.0 + 5.0 * shear) + 7.0, 4, fp * 40.0);
        float bands = 0.5 + 0.5 * sin(lat * 11.0 + 2.5 * shear);
        float glowT = 0.1 + 0.2 * bands + 1.25 * streaks * streaks + 0.12 * wisps + 8.0 * (tFactor - 1.0);
        glowT *= 0.55 + 0.45 * sqrt(mu);
        const int N = 5;
        const float stops[N] = float[N](0.0, 0.25, 0.5, 0.75, 1.0);
        const vec3 cols[N] = vec3[N](vec3(0.06, 0.01, 0.02), vec3(0.3, 0.06, 0.09), vec3(0.62, 0.17, 0.12),
                                     vec3(1.0, 0.45, 0.2), vec3(1.0, 0.72, 0.55));
        float t = clamp(glowT, 0.0, 1.0);
        vec3 c = cols[N - 1];
        for (int k = 1; k < N; ++k) {
            if (t <= stops[k]) {
                c = mix(cols[k - 1], cols[k], (t - stops[k - 1]) / (stops[k] - stops[k - 1]));
                break;
            }
        }
        color = c * (1.15 * uSurfaceBrightness);
    }
#elif SURFACE_LOOK == 3 || SURFACE_LOOK == 4
    if (uFilmColor == 1) {
        // As films and artists' impressions show them: a white-hot core fading to electric
        // blue at the limb, the surface seething with plasma (fixed to the star, so it shows
        // the rotation). The physics' brightness pattern (hot poles, hot spots) rides on it.
        float plasma = fbmFiltered(n * 5.0 + 1.7, 5, fp * 5.0) + 0.7 * ridgedFiltered(n * 9.0 + 3.0, 4, fp * 9.0);
        float I = pow(tFactor, 4.0) * (0.8 + 0.5 * plasma);
        float level = kElectricLevel;
        bool beacon = false;
#if SURFACE_LOOK == 3
        // A white dwarf: the turbulence shows across a bluer disk.
        float t = 0.42 + 0.36 * sqrt(mu) * I + 0.3 * (I - 1.0);
#else
        float t = 0.5 + 0.42 * sqrt(mu) * I + 0.12 * (I - 1.0);
        if (uBeacon == 1) {
            // As films show a spinning neutron star: a sphere of blue plasma crazed with
            // glowing veins (fixed to the star, so they show its spin), brightening toward
            // a white-blue rim where the line of sight grazes its glowing surface layer.
            vec3 cells = worley(n * 7.0, 0.0);
            float vein = 1.0 - smoothstep(0.0, 0.08 + 14.0 * fp, cells.y - cells.x);
            float fine = fbmFiltered(n * 18.0 + 4.1, 3, fp * 18.0);
            float rim = pow(1.0 - mu, 2.5);
            t = 0.12 + 0.12 * plasma + 0.45 * vein * (0.75 + 0.25 * fine) + 0.7 * rim + 0.2 * (I - 1.0);
            level = 1.15;
            beacon = true;
        }
#endif
        color = (beacon ? beaconHue(t) : electricHue(t)) * (level * uSurfaceBrightness);
    }
#endif

    vec3 hot = vec3(0.0);  // extra light from hot spots
    if (uSurfaceDetail == 1) {
        // Magnetic polar caps, heated by gas channelled onto them or by currents in the
        // magnetosphere. Their size: where the field lines that leave through the light
        // cylinder or reach the disk's inner edge land (src/app.cpp).
        if (uHotSpotTemp > 0.0) {
            float cap = max(dot(n, uMagneticAxis), dot(n, -uMagneticAxis));
            float edge = max(fp * 2.0, 0.08 * (1.0 - uCapCos));
            float spot = smoothstep(uCapCos - edge, uCapCos + edge, cap);
            if (spot > 0.0) {
                float beam = 0.0;
                if (uPulsar == 1) {
                    // The beam points along the magnetic axis: a flash when it sweeps the camera.
#if SURFACE_LOOK == 4
                    vec3 axis = magneticAxisShaken(turned + uSurfaceShift + atan(spinA, R), delay);
#else
                    vec3 axis = magneticAxisAt(turned + uSurfaceShift + atan(spinA, R));
#endif
                    beam = exp(-(1.0 - abs(dot(-rayDir, axis))) / 0.01);
                }
                float Tcap = uHotSpotTemp * pow(limbTemperature(mu), tExp);
                hot += spot * shiftedBlackbody(Tcap, g, Tref) * kPhotosphereLevel * uSurfaceBrightness * (1.0 + 40.0 * beam);
            }
        }
        // Bodies that struck the surface, glowing where they hit and fading.
        for (int k = 0; k < uImpactCount; ++k) {
            float age = uImpactAges[k] - delay;
            if (age < 0.0) continue;  // its light has not reached us yet
            float strength = uImpacts[k].w * exp(-age / uImpactFade);
            float spot = exp(-(1.0 - dot(n, uImpacts[k].xyz)) * 300.0);
            hot += strength * spot * shiftedBlackbody(max(uSurfaceTemp * 2.5, 20000.0), g, uSurfaceTemp);
        }
#if SURFACE_LOOK == 1
        // White-light flares: a small patch of the surface flashes at ~9,500 K, rising
        // in about a minute and fading over ten (Kowalski et al. 2013).
        for (int k = 0; k < uFlareCount; ++k) {
            float age = uFlareAges[k] - delay;
            if (age < 0.0) continue;  // its light has not reached us yet
            float curve = age < uFlareRise ? age / uFlareRise : exp(-(age - uFlareRise) / uFlareDecay);
            float d2 = max(2.0 - 2.0 * dot(n, uFlares[k].xyz), 0.0);
            float kernel = exp(-d2 / (0.0012 + fp * fp));
            hot += uFlares[k].w * curve * kernel * shiftedBlackbody(9500.0, g, Tref) * kPhotosphereLevel * uSurfaceBrightness;
        }
#endif
        // Gas that struck it (a map in the rotating frame, drawn as it is now).
        float gas = texture(uImpactMap, vec2(th / PI, fract(lon / (2.0 * PI)))).r;
        if (gas > 1e-4) hot += min(gas, 50.0) * shiftedBlackbody(max(uSurfaceTemp * 2.5, 20000.0), g, uSurfaceTemp);
    }
    color += hot * (2.0 + 3.0 * mu) / 5.0;
    return color * (1.0 + flash);
}

#if SURFACE_LOOK == 4
// --- A neutron star's magnetosphere ---------------------------------------------------
//
// A dipole: field lines r = L sin^2(theta_m) about the magnetic axis (theta_m the
// magnetic colatitude, L a line's apex). Lines whose apex lies inside the light cylinder
// c / Omega (or inside the disk's inner edge, where the disk's gas takes over) close back
// on the star and turn with it; plasma held on them glows (a magnetar's, carrying the
// currents of its twisted field, in X-rays). Those beyond are dragged open and stream
// pair plasma out over the poles, where a pulsar's beams are made. Drawn in false color
// (X-ray and radio light shown as blue): none of it is bright in visible light.
uniform float uFieldGlow;     // brightness of the closed field lines (0: not drawn)
uniform float uFieldTwist;    // a magnetar's twist: azimuth between a line's footpoints, radians
uniform float uFieldReach;    // field lines are drawn out to this apex (render units)
uniform float uPairGlow;      // a magnetar's pair corona, the diffuse glow of its twisted field
uniform float uBeamGlow;      // brightness of the beams (0: none)
uniform float uBeamWidth;     // their half-opening angle, radians

const float kFieldFirst = 1.6;   // the innermost drawn line's apex, in surface radii
const float kFieldRatio = 1.45;  // apex ratio between successive shells of lines
const float kFieldLines = 9.0;   // lines per shell around the axis
const float kFieldWidth = 0.022; // a line's half-width, in surface radii
const vec3 kFieldColor = vec3(0.42, 0.62, 1.0);
const vec3 kTwistedColor = vec3(0.62, 0.5, 1.0);
const vec3 kBeamColor = vec3(0.62, 0.8, 1.0);
const float kHaloLevel = 0.12;  // brightness of a young neutron star's glow, with uHaloGlow 1
const float kRimLevel = 1.6;   // ... and of the thin glowing layer that rings it

// A point's place among the field lines: returns the azimuth of its line in units of
// the spacing between lines, and L and cos(theta_m). A twisted line's azimuth runs from
// -twist/2 at one footpoint to +twist/2 at the other; this undoes that.
float fieldCoordinate(vec3 x, vec3 axis, vec3 e1, vec3 e2, float lines, out float L, out float cm) {
    float r = length(x);
    cm = dot(x, axis) / r;
    L = r / max(1.0 - cm * cm, 1e-6);
    float phi = atan(dot(x, e2), dot(x, e1));
    float cFoot = sqrt(max(1.0 - uSurfaceR / L, 1e-3));
    phi -= 0.5 * uFieldTwist * clamp(cm / cFoot, -1.0, 1.0);
    return phi * lines / (2.0 * PI);
}

// The field lines, pair corona and beams along the chord X0 -> X1, which the light left
// delay0 and delay1 ago; pathLength: how far it then had to go to the camera.
void traceMagnetosphere(vec3 X0, vec3 X1, float delay0, float delay1, float pathLength, inout vec3 color,
                        float transmittance) {
    float R = uSurfaceR;
    vec3 d = X1 - X0;
    float chord = length(d);
    if (chord <= 0.0) return;
    vec3 dir = d / chord;
    // The chord's closest approach to the star.
    float tc = clamp(-dot(X0, dir), 0.0, chord);
    float rClosest = length(X0 + tc * dir);

    // The field turns rigidly with the star (inside the light cylinder its pattern is
    // quasi-static): its axis when the light left the chord's middle.
    float delayMid = 0.5 * (delay0 + delay1);
    float turned = uSurfacePhase - uSurfaceOmega * delayMid;
    vec3 axis = magneticAxisAt(turned + uSurfaceShift + atan(spinA, rClosest));
    // The magnetic frame; pinned to the star's turning even for an aligned axis.
    vec3 e1 = abs(axis.z) < 0.999 ? normalize(cross(vec3(0.0, 0.0, 1.0), axis)) : vec3(cos(turned), sin(turned), 0.0);
    vec3 e2 = cross(axis, e1);

    if (uFieldGlow > 0.0 && rClosest < uFieldReach) {
        // Lines are sheets of constant (twisted) magnetic azimuth crossed by the chord:
        // find each crossing (half-integer coordinates: alternate shells are staggered by
        // half a spacing) and measure, within the sheet, how far it passes from the line.
        // A Gaussian tube of half-width w seen across collects sqrt(pi) w of its glow, more
        // when the ray runs along its sheet. No sampling: as sharp for any step length.
        // In film colors, as artists draw them: denser shells of more lines, each with a
        // soft halo, gently wavering.
        bool film = uFilmColor == 1;
        float lines = film ? (uFieldTwist > 0.0 ? 16.0 : 12.0) : kFieldLines;
        float ratio = film ? 1.3 : kFieldRatio;
        float L0, cm0, L1, cm1;
        float c0 = fieldCoordinate(X0, axis, e1, e2, lines, L0, cm0);
        float c1 = fieldCoordinate(X1, axis, e1, e2, lines, L1, cm1);
        float dc = c1 - c0;
        dc -= lines * floor(dc / lines + 0.5);
        if (abs(dc) > 1e-6) {
            float stepSign = dc > 0.0 ? 1.0 : -1.0;
            float first = dc > 0.0 ? floor(2.0 * c0) + 1.0 : ceil(2.0 * c0) - 1.0;  // in half spacings
            for (int k = 0; k < (film ? 8 : 6); ++k) {
                float target = 0.5 * (first + stepSign * float(k));
                float t = (target - c0) / dc;
                if (t < 0.0 || t > 1.0) break;
                vec3 x = X0 + t * d;
                float r = length(x);
                if (r < R) continue;
                float L, cm;
                fieldCoordinate(x, axis, e1, e2, lines, L, cm);
                float u = log(L / (kFieldFirst * R)) / log(ratio);
                float shell = floor(u + 0.5);
                if (shell < 0.0) continue;
                float Lk = kFieldFirst * R * pow(ratio, shell);
                float lineIndex = floor(mod(target + 1000.0 * lines, lines) * 2.0 + 0.5);
                if (film) Lk *= 1.0 + 0.035 * gradientNoise(vec3(cm * 4.0, lineIndex * 1.37, shell * 2.1));
                if (Lk > uFieldReach) continue;
                // Alternate shells hold the integer and the half-integer sheets.
                float halfStep = mod(2.0 * target, 2.0);
                if (abs(halfStep - mod(shell, 2.0)) > 0.5) continue;
                float s2 = max(1.0 - cm * cm, 1e-6);
                // Within a meridional plane, neighbouring lines lie dL sin^3 / sqrt(1 + 3 cos^2) apart.
                float gap = abs(L - Lk) * s2 * sqrt(s2) * inversesqrt(1.0 + 3.0 * cm * cm);
                float w0 = kFieldWidth * R * sqrt(Lk / R);
                float w = max(w0, 0.7 * uPixelAngle * (pathLength + t * chord));
                float profile = exp(-gap * gap / (w * w)) * (w0 / w);
                if (film) profile += 0.2 * exp(-gap * gap / (9.0 * w * w)) * (w0 / w);
                if (profile < 1e-3) continue;
                // A ray running along the sheet stays in the line longer.
                vec3 eAz = normalize(cross(axis, x));
                float along = 1.0 / max(abs(dot(dir, eAz)), 0.2);
                // Each line its own brightness; brighter near the star, where the field is.
                float vary = 0.45 + 0.9 * hash13(vec3(lineIndex, shell, 7.3));
                float strength = uFieldGlow * vary * pow(R / r, 1.5) * profile * along * 1.7725 * w0 / R;
                vec3 tint = mix(kFieldColor, kTwistedColor, clamp(uFieldTwist, 0.0, 1.0));
                // Electric blue, white-hot near the star.
                if (film) tint = 2.0 * mix(vec3(0.12, 0.42, 1.0), vec3(0.75, 0.95, 1.0), clamp(1.6 * R / r - 0.4, 0.0, 1.0));
                color += transmittance * tint * strength;
            }
        }
    }

    if (uPairGlow > 0.0 && rClosest < 6.0 * R) {
        // A magnetar's twisted field is threaded by currents of electrons and positrons,
        // which scatter the surface's X-rays: a glow falling off steeply (~ r^-4).
        const int N = 3;
        for (int k = 0; k < N; ++k) {
            float f = (float(k) + 0.5) / float(N);
            vec3 x = X0 + f * d;
            float r = length(x);
            if (r < R) continue;
            float cm = dot(x, axis) / r;
            // Densest on the twisted lines over the magnetic equator and mid-latitudes.
            float shape = 0.6 + 0.4 * (1.0 - cm * cm);
            float q = R / r;
            vec3 tint = uFilmColor == 1 ? vec3(0.2, 0.5, 1.6) : kTwistedColor;
            color += transmittance * tint * (uPairGlow * shape * q * q * q * q * chord / float(N) / R);
        }
    }

    if (uHaloGlow > 0.0 && rClosest < 12.0 * R) {
        // The glow round a young neutron star: its X-rays scattered by the plasma held
        // in its magnetosphere and streaming off its poles, falling off as r^-3, white
        // near the star and blue farther out; a glitch's flash brightens it.
        float since = uFlash.x - delayMid;
        float flare = since >= 0.0 ? uFlash.y * exp(-since / uFlash.z) : 0.0;
        const int N = 4;
        for (int k = 0; k < N; ++k) {
            vec3 x = X0 + ((float(k) + 0.5) / float(N)) * d;
            float r = length(x);
            if (r < R) continue;
            float q = R / r;
            vec3 tint = mix(vec3(0.03, 0.3, 1.0), vec3(0.15, 0.7, 1.0), q * q);
            // A thin glowing layer at the surface: seen edge-on it rings the star.
            float shell = kRimLevel * exp(-(r - R) / (0.08 * R));
            color += transmittance * (tint * kHaloLevel * q * q * q + vec3(0.2, 0.75, 1.0) * shell) *
                     (uHaloGlow * (1.0 + 2.0 * flare) * chord / float(N) / R);
        }
    }

    if (uBeamGlow > 0.0) {
        // The beams, shown as if they lit a faint haze around the star (a searchlight in
        // fog): flux ~ 1 / r^2 inside a cone of half-angle uBeamWidth about each magnetic
        // pole. Light seen at distance r left the star r / c earlier, when the star had
        // turned less: the beams curve back into a spiral (strongly so for millisecond
        // pulsars, whose light cylinder is a few radii out). In reality a beam is seen
        // only as the flash when it sweeps across you.
        float rMid = max(length(X0 + 0.5 * d), R);
        // Skip chords that pass nowhere near either cone: seen from the star a chord spans
        // the angle between its ends, and every point of it lies within half that of an
        // end or the middle; the beam turns by Omega times the extra delay along it.
        vec3 axisMid = magneticAxisAt(uSurfacePhase - uSurfaceOmega * (delayMid + rMid - R) + uSurfaceShift + atan(spinA, rMid));
        vec3 u0 = normalize(X0), u1 = normalize(X1), um = normalize(X0 + 0.5 * d);
        float span = acos(clamp(dot(u0, u1), -1.0, 1.0));
        float closest = max(abs(dot(u0, axisMid)), max(abs(dot(u1, axisMid)), abs(dot(um, axisMid))));
        // A beacon's needle sits in a glow; starquakes swing the poles.
        float reach = (uBeacon == 1 ? 5.3 : 1.3) * uBeamWidth + 0.5 * span + abs(uSurfaceOmega) * chord +
                      (uQuakeCount > 0 ? 0.12 : 0.0);
        if (acos(min(closest, 1.0)) > reach) return;  // the beams come last
        int n = clamp(int(ceil(chord / (0.3 * rMid * uBeamWidth))), 1, 8);
        vec3 cameraX = toCartesian(uCamPos.x, uCamPos.y, uCamPos.z);
        float dl = chord / float(n);
        for (int k = 0; k < n; ++k) {
            float f = (float(k) + 0.5) / float(n);
            vec3 x = X0 + f * d;
            float r = length(x);
            if (r < R) continue;
            float delay = mix(delay0, delay1, f) + (r - R);
            vec3 beamAxis = magneticAxisShaken(uSurfacePhase - uSurfaceOmega * delay + uSurfaceShift + atan(spinA, r), delay);
            float along = dot(x, beamAxis);
            float c = abs(along) / r;
            float psi = acos(min(c, 1.0)) / uBeamWidth;
            float profile, glowPart = 0.0;
            if (uBeacon == 1) {
                // The core beam, a needle of half width uBeamWidth at half maximum, in a soft
                // glow. Each pulse is brighter or fainter than the last, and flickers within
                // (as real pulsars' do; their radio light comes in bursts down to
                // microseconds long). Set by when the light left the star, so each change
                // runs out along the beam at the speed of light.
                // ... and a white-hot thread along its middle.
                profile = exp(-0.693 * psi * psi) + 0.7 * exp(-12.0 * psi * psi);
                glowPart = 0.06 * exp(-0.173 * psi * psi);
                if (profile + glowPart < 1e-3) continue;
                float turns = uSurfaceTurns - uSurfaceOmega * delay / (2.0 * PI);
                float side = along > 0.0 ? 0.0 : 17.0;
                float k0 = floor(turns);
                float pulse = mix(hash13(vec3(k0, side, 3.1)), hash13(vec3(k0 + 1.0, side, 3.1)),
                                  smoothstep(0.0, 1.0, turns - k0));
                float flicker = 0.85 + 0.15 * sin(2.0 * PI * 3.1 * turns + side) * sin(2.0 * PI * 4.7 * turns + 1.3);
                float strength = (0.45 + 1.1 * pulse) * flicker;
                profile *= strength;
                glowPart *= strength;
            } else {
                // A core inside a hollow cone, as in most radio pulsars: the core's half width
                // at half maximum is 1.22 deg / sqrt(P), a fifth of the cone's 5.75 deg / sqrt(P)
                // (Rankin 1990, 1993); the cone's ring lies just inside the beam's edge.
                float cone = (psi - 0.85) / 0.12;
                profile = exp(-psi * psi * 15.7) + 0.35 * exp(-cone * cone);
                if (profile < 1e-3) continue;
            }
            float q = R / r;
            float rise = smoothstep(R, 1.6 * R, r);
            vec3 tint = kBeamColor;
            // Films' beams: white-hot along the axis, electric blue at the edges.
            if (uFilmColor == 1) tint = 2.0 * mix(vec3(0.25, 0.6, 1.0), vec3(0.95, 0.98, 1.0), exp(-psi * psi * 6.0));
            // A beacon's needle: cyan-white, its edges deep blue.
            if (uFilmColor == 1 && uBeacon == 1) {
                tint = 2.0 * mix(vec3(0.05, 0.45, 1.0), vec3(0.55, 0.95, 1.0), exp(-psi * psi * 2.0));
                tint = mix(tint, vec3(2.2), exp(-psi * psi * 14.0));
            }
            // Films keep the needle bright far out: it fades as 1 / r, not 1 / r^2 (its
            // soft glow as 1 / r^2), out to 30 radii; it fades near the camera, so a beam
            // turned toward you does not swell into a blob as it passes the lens.
            float dim = uBeacon == 1 && uFilmColor == 1 ? q : q * q;
            if (uBeacon == 1) {
                dim *= smoothstep(0.15 * uCamPos.x, 0.4 * uCamPos.x, length(x - cameraX)) *
                       (1.0 - smoothstep(25.0 * R, 30.0 * R, r));
            }
            // Its soft glow with a violet fringe, as in films.
            vec3 glowTint = uFilmColor == 1 ? vec3(0.9, 0.8, 2.0) : tint;
            color += transmittance * (tint * profile * dim + glowTint * glowPart * q * q * smoothstep(0.0, 1.0, dim / max(q, 1e-6))) *
                     (uBeamGlow * rise * dl / R);
        }
    }
}
#endif

#if SURFACE_LOOK == 3
// --- A cataclysmic variable's companion and stream ----------------------------------------
//
// A red dwarf filling its Roche lobe (Eggleton's radius; drawn as an ellipsoid stretched
// toward L1), and the stream of its gas falling from L1 to the disk's rim along the
// ballistic path of the Roche potential (src/central.cpp, rocheStream), bent ahead of
// the companion by the Coriolis force. Both turn with the binary. Light bends by ~1e-5
// rad around a white dwarf: their shapes are met along the (straight) ray chords.
uniform vec4 uCompanion;       // center (render units) and Roche-lobe radius (0: none)
uniform vec3 uCompanionAxis;   // from the white dwarf toward it
uniform float uCompanionTemp;  // its surface temperature, K
uniform int uStreamCount;
uniform vec3 uStream[16];      // the stream from L1 to the rim (render units)
uniform float uStreamWidth;    // its half-width at L1 (render units); it widens toward the disk

// Semi-axes over the lobe radius (toward L1, along the orbit, up), centered 0.1 of it
// toward L1: its tip reaches L1, 1.35 lobe radii out for q = 0.25.
const vec3 kLobeShape = vec3(1.25, 0.95, 0.93);
const float kLobeShift = 0.1;

// The companion's frame: toward the white dwarf's far side (x), along the orbit, up.
void companionFrame(out vec3 ax, out vec3 ay, out vec3 az) {
    ax = uCompanionAxis;
    az = vec3(0.0, 0.0, 1.0);
    ay = normalize(cross(az, ax));
}

// Where the chord X0 -> X1 enters the companion, as a fraction along it (2 if not).
float companionCrossing(vec3 X0, vec3 X1) {
    vec3 ax, ay, az;
    companionFrame(ax, ay, az);
    vec3 c = uCompanion.xyz - kLobeShift * uCompanion.w * ax;  // the lobe bulges toward L1
    vec3 k = 1.0 / (kLobeShape * uCompanion.w);
    vec3 p = X0 - c, dd = X1 - X0;
    vec3 a = vec3(dot(p, ax), dot(p, ay), dot(p, az)) * k;
    vec3 b = vec3(dot(dd, ax), dot(dd, ay), dot(dd, az)) * k;
    float bb = dot(b, b), ab = dot(a, b), cc = dot(a, a) - 1.0;
    float disc = ab * ab - bb * cc;
    if (bb <= 0.0 || disc < 0.0) return 2.0;
    float t = (-ab - sqrt(disc)) / bb;
    return t >= 0.0 && t <= 1.0 ? t : 2.0;
}

// The companion's light at its surface point p, seen along rayDir. A red dwarf, darker
// at its limb, its convection mottling it, and the hemisphere facing the white dwarf and
// its disk heated by their light (irradiation: ~20% hotter there).
vec3 companionRadiance(vec3 p, vec3 rayDir, float dist) {
    vec3 ax, ay, az;
    companionFrame(ax, ay, az);
    vec3 c = uCompanion.xyz - kLobeShift * uCompanion.w * ax;
    vec3 q = p - c;
    vec3 local = vec3(dot(q, ax), dot(q, ay), dot(q, az)) / (kLobeShape * uCompanion.w);
    vec3 normal = normalize(ax * (local.x / kLobeShape.x) + ay * (local.y / kLobeShape.y) + az * (local.z / kLobeShape.z));
    float mu = clamp(dot(normal, -rayDir), 0.0, 1.0);
    float facing = max(dot(normal, -ax), 0.0);  // toward the white dwarf
    float fp = dist * uPixelAngle / uCompanion.w;
    vec3 nl = normalize(local);  // tidally locked: the pattern is fixed in its frame
    float cells = fbmFiltered(nl * 9.0 + 4.0, 5, fp * 9.0) + 0.6 * ridgedFiltered(nl * 18.0, 4, fp * 18.0);
    float T = uCompanionTemp * (1.0 + 0.2 * facing * facing) * (0.97 + 0.04 * cells) * limbTemperature(mu);
    if (uFilmColor == 1) {
        // Films' red dwarf, in the Sun's ultraviolet colors made redder: boiling, bright
        // toward its limb, its face toward the white dwarf heated orange-yellow, and
        // brightest where its gas pours out at L1.
        float flames = ridgedFiltered(nl * 6.0 + 9.0, 5, fp * 6.0);
        float t = 0.3 + 0.16 * cells + 0.3 * flames * flames + 0.18 * pow(1.0 - mu, 3.0) + 0.22 * facing * facing +
                  0.35 * pow(max(dot(nl, -vec3(1.0, 0.0, 0.0)), 0.0), 12.0);
        return aia304(t) * (1.25 * uSurfaceBrightness);
    }
    return shiftedBlackbody(T, 1.0, uSurfaceTemp * kCenterTemperature) * (kPhotosphereLevel * uSurfaceBrightness);
}

// The stream's glow and shadow along the chord X0 -> X1: a tube around the stream's
// polyline, half-width uStreamWidth at L1 growing threefold to the rim, its gas heated
// from ~4,000 K at L1 to ~12,000 K as it falls (the shock at the rim is the disk's
// bright spot). Locally straight: sqrt(pi) w / sin(gamma) exp(-D^2 / w^2) across it.
void traceStream(vec3 X0, vec3 X1, float pathLength, inout vec3 color, inout float transmittance) {
    if (uStreamCount < 2) return;
    vec3 d1 = X1 - X0;
    float a11 = dot(d1, d1);
    float chord = sqrt(a11);
    if (chord <= 0.0) return;
    // Quick reject: the chord passes far from the stream's bounding sphere.
    vec3 mid = 0.5 * (uStream[0] + uStream[uStreamCount - 1]);
    float bound = 0.5 * length(uStream[0] - uStream[uStreamCount - 1]) + 4.0 * uStreamWidth + 0.3 * length(uStream[0]);
    float tc = clamp(dot(mid - X0, d1) / a11, 0.0, 1.0);
    if (length(X0 + tc * d1 - mid) > bound) return;
    // The segment the chord passes closest to.
    float bestD2 = 1e30, bestF = 0.0, bestU = 0.0, bestSin = 1.0;
    for (int i = 0; i + 1 < uStreamCount; ++i) {
        vec3 sa = uStream[i];
        vec3 d2 = uStream[i + 1] - sa;
        vec3 w0v = X0 - sa;
        float a22 = dot(d2, d2), a12 = dot(d1, d2);
        float b1 = dot(d1, w0v), b2 = dot(d2, w0v);
        float den = a11 * a22 - a12 * a12;
        float u = den > 1e-12 ? clamp((a12 * b2 - a22 * b1) / den, 0.0, 1.0) : 0.0;
        float v = a22 > 1e-12 ? clamp((a12 * u + b2) / a22, 0.0, 1.0) : 0.0;
        u = clamp((a12 * v - b1) / a11, 0.0, 1.0);
        vec3 gap = (X0 + u * d1) - (sa + v * d2);
        float D2 = dot(gap, gap);
        if (D2 < bestD2) {
            bestD2 = D2;
            bestF = (float(i) + v) / float(uStreamCount - 1);  // 0 at L1, 1 at the rim
            bestU = u;
            bestSin = max(length(cross(d1 / chord, d2 * inversesqrt(max(a22, 1e-12)))), 0.2);
        }
    }
    float f = bestF;
    // It leaves L1 thin and spreads as it falls.
    float w0 = uStreamWidth * (1.0 + 2.5 * f * f);
    float w = max(w0, 0.7 * uPixelAngle * (pathLength + bestU * chord));
    if (bestD2 > 9.0 * w * w) return;
    // Knots and wisps along it.
    float clumps = 0.6 + 0.6 * smoothstep(-0.4, 0.5, gradientNoise(vec3(f * 18.0, 3.1, 0.0)));
    float depth = 1.7725 * exp(-bestD2 / (w * w)) / bestSin * (w0 / w);
    float alpha = 1.0 - exp(-depth * clumps);
    float T = mix(3500.0, 14000.0, f * f);
    vec3 glow = uFilmColor == 1 ? aia304(0.42 + 0.55 * f) * (0.9 + 0.6 * f)
                                : shiftedBlackbody(T, 1.0, uSurfaceTemp * kCenterTemperature) * kPhotosphereLevel;
    color += transmittance * glow * (alpha * uSurfaceBrightness);
    transmittance *= 1.0 - 0.8 * alpha;
}
#endif

#if SURFACE_LOOK == 4
// --- Helical fields around a neutron star's jets --------------------------------------------
//
// The disk's rotation winds the field lines it launches into helices around each jet (as
// in simulations of magnetized jets, and as their polarization shows). Three strands per
// jet on its surface, their winding loosening outward. Drawn exactly where the chord
// crosses each strand's sheet (constant azimuth less the winding), like the field lines.
const float kJetStrands = 3.0;

float helixCoordinate(vec3 x) {
    float z = abs(x.z);
    float winding = 7.0 * log(1.0 + z / 6.0);  // radians: a turn every ~6 M near the base
    return (atan(x.y, x.x) - sign(x.z) * winding) * kJetStrands / (2.0 * PI);
}

void traceJetHelices(vec3 X0, vec3 X1, float pathLength, inout vec3 color, float transmittance) {
    const float zMax = 70.0;
    float z0 = uSurfaceR + 1.0;
    if (min(abs(X0.z), abs(X1.z)) > zMax || max(abs(X0.z), abs(X1.z)) < z0) return;
    // Skip chords far from the axis.
    vec2 a = X0.xy, dxy = X1.xy - X0.xy;
    float tc = clamp(-dot(a, dxy) / max(dot(dxy, dxy), 1e-12), 0.0, 1.0);
    float zm = max(abs(X0.z), abs(X1.z));
    if (length(a + tc * dxy) > 2.0 * (0.4 + zm * uJetOpening) + 1.0) return;
    vec3 d = X1 - X0;
    float chord = length(d);
    float c0 = helixCoordinate(X0), c1 = helixCoordinate(X1);
    float dc = c1 - c0;
    dc -= kJetStrands * floor(dc / kJetStrands + 0.5);
    if (abs(dc) < 1e-6) return;
    float stepSign = dc > 0.0 ? 1.0 : -1.0;
    float first = dc > 0.0 ? floor(c0) + 1.0 : ceil(c0) - 1.0;
    for (int k = 0; k < 6; ++k) {
        float t = (first + stepSign * float(k) - c0) / dc;
        if (t < 0.0 || t > 1.0) break;
        vec3 x = X0 + t * d;
        float z = abs(x.z);
        if (z < z0 || z > zMax) continue;
        float w = 0.4 + z * uJetOpening;
        float gap = length(x.xy) - 0.85 * w;
        float lw0 = 0.05 + 0.01 * z;
        float lw = max(lw0, 0.7 * uPixelAngle * (pathLength + t * chord));
        float profile = (exp(-gap * gap / (lw * lw)) + 0.25 * exp(-gap * gap / (9.0 * lw * lw))) * (lw0 / lw);
        if (profile < 1e-3) continue;
        float fade = pow(z0 / z, 0.6) * smoothstep(z0, z0 + 2.0, z) * (1.0 - smoothstep(0.7 * zMax, zMax, z));
        color += transmittance * vec3(0.3, 0.9, 1.0) * (2.2 * uJetPower * profile * fade);
    }
}
#endif

#if SURFACE_LOOK == 4
// --- The thick inner flow around an accreting neutron star ----------------------------------
//
// Near a bright neutron star the disk is no thin sheet: radiation pressure and heating
// puff its inner part up (H / r ~ 0.1-0.3 near the Eddington rate), and simulations show
// a thick, turbulent torus of plasma whose eddies the Keplerian shear winds into spiral
// flames, hottest inward. Drawn (film colors) as a glowing, absorbing volume along the
// bent rays, its flow seen where it was when the light left it, its approaching side
// beamed brighter (g^3, optically thin light).
uniform float uFlowGlow;  // its brightness (0: none)

// Three octaves of value noise: the flow's eddies (cheap: it is read many times per ray).
float flowNoise(vec3 p) {
    return (0.5 * valueNoise(p) + 0.25 * valueNoise(p * 2.07 + 13.7) + 0.125 * valueNoise(p * 4.3 + 7.3)) / 0.875;
}

float flowDensity(vec3 x, float r, float time) {
    float H = 0.17 * r;
    float zz = x.z / H;
    float vertical = exp(-0.5 * zz * zz);
    if (vertical < 0.02) return 0.0;
    float edges = smoothstep(uDiskInner * 0.9, uDiskInner * 1.3, r) * (1.0 - smoothstep(0.65 * uDiskOuter, uDiskOuter, r));
    if (edges <= 0.0) return 0.0;
    // Two copies of the eddies, each wound by the shear over one flow period, cross-faded
    // (as for the thin disk) so the spirals renew before they wind too tight.
    float omega = uDiskDirection * sqrt(MASS) / (r * sqrt(r));
    float phi = atan(x.y, x.x);
    float phase1 = fract(time / uFlowPeriod), phase2 = fract(phase1 + 0.5);
    float lr = log(r);
    // Trailing spiral streams: the eddies' azimuth runs on with ln r.
    float a1 = phi - omega * phase1 * uFlowPeriod + 2.4 * lr, a2 = phi - omega * phase2 * uFlowPeriod + 2.4 * lr;
    vec3 q1 = vec3(cos(a1) * 3.5, sin(a1) * 3.5, lr * 2.6 + zz * 0.6);
    vec3 q2 = vec3(cos(a2) * 3.5, sin(a2) * 3.5, lr * 2.6 + zz * 0.6 + 17.0);
    float n = mix(flowNoise(q1 * 1.8), flowNoise(q2 * 1.8), abs(2.0 * phase1 - 1.0));
    // Wispy at the surface, denser in the core, broken into streaks.
    float flames = smoothstep(0.42 + 0.18 * abs(zz), 0.85, n);
    return vertical * edges * flames;
}

void traceThickFlow(vec3 X0, vec3 X1, float delay0, float delay1, inout vec3 color, inout float transmittance) {
    // Quick rejects: beyond the disk, or well above or below it.
    float rc0 = length(X0.xy), rc1 = length(X1.xy);
    float rMax = max(rc0, rc1);
    if (min(rc0, rc1) > uDiskOuter * 1.05 && min(length(X0), length(X1)) > uDiskOuter * 1.05) return;
    float zLimit = 0.6 * max(rMax, uDiskInner);
    if (X0.z * X1.z > 0.0 && min(abs(X0.z), abs(X1.z)) > zLimit) return;
    float chord = length(X1 - X0);
    float rNear = max(min(length(X0), length(X1)), uDiskInner);
    int n = clamp(int(ceil(chord / (0.15 * rNear))), 2, 4);
    float dl = chord / float(n);
    float jitter = fract(52.9829189 * fract(dot(gl_FragCoord.xy, vec2(0.06711056, 0.00583715))));
    // The flow's frequency shift, once per chord (it varies slowly along it): circular
    // orbits at the chord's middle radius (Kerr, prograde).
    float rMid = max(length(mix(X0, X1, 0.5).xy), uDiskInner);
    float sm = sqrt(MASS);
    float omega = uDiskDirection * sm / (pow(rMid, 1.5) + uDiskDirection * spinA * sm);
    float gtt, gtp, gpp;
    blMetric(rMid, 0.5 * PI, gtt, gtp, gpp);
    float norm = -(gtt + 2.0 * gtp * omega + gpp * omega * omega);
    float g = norm > 0.0 ? 1.0 / (inversesqrt(norm) * (photonE - photonL * omega)) : 1.0;
    for (int k = 0; k < n; ++k) {
        float f = (float(k) + jitter) / float(n);
        vec3 x = mix(X0, X1, f);
        float r = length(x.xy);
        if (length(x) < uSurfaceR) continue;
        float delay = mix(delay0, delay1, f);
        float density = flowDensity(x, r, uFlowTime - delay);
        if (density < 1e-3) continue;
        float heat = pow(uDiskInner / r, 0.75);
        float t = (0.22 + 0.7 * heat) * (0.7 + 0.5 * density) * sqrt(clamp(g, 0.3, 3.0));
        float alpha = 1.0 - exp(-2.6 * density * dl / r);
        color += transmittance * aia304(t) * (uFlowGlow * alpha * clamp(g * g * g, 0.05, 8.0) * (0.6 + 2.4 * heat * heat));
        transmittance *= 1.0 - alpha;
        if (transmittance < 0.01) return;
    }
}
#endif
