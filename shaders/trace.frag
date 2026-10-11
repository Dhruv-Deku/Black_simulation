#version 460 core

// Kerr black hole ray tracer.
//
// Each pixel follows the photon that reaches the camera, traced backwards along an
// exact null geodesic of the Kerr metric in ingoing Kerr-Schild coordinates
// (src/kerr.h), which stay regular across the horizon: the camera can be inside it.
// Units are G = c = M = 1. The photon is scaled so the camera measures unit energy;
// its constants of motion are E, L and Q (E can be negative in the ergosphere).
//
// The equations are Carter's separated geodesic equations in Mino time tau
// (d lambda = Sigma d tau), used in second-order form so that turning points need
// no special handling (see src/kerr.h; src/kerr.cpp is the tested CPU reference):
//     d^2 r / d tau^2     = R'(r) / 2
//     d^2 theta / d tau^2 = Theta'(theta) / 2
//     d phi / d tau       = a F - a E + L / sin^2(theta)
//     d t / d tau         = (r^2 + a^2) F - a (a E sin^2(theta) - L) - dr/dtau
// with F = (P + dr/dtau) / Delta, evaluated without the 0/0 at the horizon.
// Because we trace backwards, s = -tau: the second derivatives keep their sign and
// the first-order rates flip. Integrating dt/dtau along the way gives the
// light-travel delay: how long ago the light was emitted at each point of the ray.

layout(location = 0) out vec4 outColor;

uniform vec2 uResolution;
uniform float uSpin;
uniform vec3 uCamPos;       // camera (r, theta, phi)
uniform mat3 uCamBasis;     // camera right/up/forward on the camera frame's (r, theta, phi) axes
// The camera frame (4-velocity and three axes), each column holding one frame vector
// as (u_t, u_phi, u_theta, u^r): a photon arriving from direction n has
// (p_t, p_phi, p_theta, p^r) = uCamPhoton * (1, -n).
uniform mat4 uCamPhoton;
uniform float uTanHalfFov;
uniform float uPixelAngle;
uniform int uMaxSteps;
uniform float uStepSize;
uniform float uEscapeRadius;

uniform int uDiskEnabled;
uniform float uDiskInner;
uniform float uDiskOuter;
uniform float uDiskTemp;
uniform float uDiskBrightness;
uniform float uDiskOpacity;
uniform float uDiskTurbulence;
uniform float uDiskDirection;  // +1 prograde, -1 retrograde
uniform float uFlowTime;       // coordinate time modulo uFlowPeriod
uniform float uFlowPeriod;
uniform int uRelativity;  // 0 = Doppler + gravitational shift, 1 = beaming only, 2 = off

uniform int uSkyMode;  // 0 = stars, 1 = grid, 2 = image, 3 = none
uniform float uStarBrightness;
uniform float uGalaxyBrightness;
uniform float uSkyImageBrightness;
uniform float uSkyLod;
uniform vec2 uBlackbodyRange;  // log10 of the first and last temperature in the table
uniform int uDebugView;

uniform int uDebrisEnabled;
uniform vec2 uGridLogR;  // log of the debris grid's inner and outer radius
uniform float uDebrisBrightness;
uniform float uDebrisOpacity;

uniform int uPathsEnabled;      // draw orbit trails and predictions as lensed lines
uniform ivec3 uPathGridDims;    // cells binning the line segments (log r, theta, phi)
uniform float uPathWidth;       // line half-width in pixels
uniform int uSkipBody;          // the body the camera rides with (not drawn), or -1

uniform int uJetsEnabled;
uniform float uJetPower;
uniform float uJetSpeed;        // beta of the jet plasma
uniform float uJetOpening;      // tan of the half opening angle
uniform float uTime;            // coordinate time, for moving knots in the jets

layout(binding = 0) uniform sampler2D uBlackbody;
layout(binding = 1) uniform sampler2D uSkyTex;
layout(binding = 2) uniform sampler2D uHistory;         // past body positions, 4-velocities and shapes
layout(binding = 3) uniform sampler3D uDebrisDensity;   // rho, rho T / 1e4, rho u^t, rho u^r
layout(binding = 4) uniform sampler3D uDebrisVelocity;  // rho u^theta, rho u^phi
layout(binding = 5) uniform usampler3D uDebrisOccupancy;  // nonzero where any debris is nearby
layout(binding = 6) uniform usampler3D uPathOccupancy;    // nonzero for path cells holding any segment

struct Body {
    vec4 positionRadius;  // Cartesian center now, radius
    vec4 fourVelocity;    // u^t, u^r, u^theta, u^phi now
    vec4 params;          // temperature, brightness, kind, seed
    vec4 motion;          // rotation angle now, rotation rate, history slot, history length
    vec4 range;           // r range over the light-delay window, time since lost (-1: it is not), bounding radius
    vec4 shapeA;          // S_xx, S_yy, S_zz, stretch: S maps the (tidally deformed) body onto a sphere
    vec4 shapeB;          // S_xy, S_xz, S_yz, fastest it moves (coordinate speed) over the window
};
const int kHistoryRows = 4;  // per body: position, 4-velocity, shapeA, shapeB
layout(std140, binding = 0) uniform Bodies {
    Body bodies[32];
    ivec4 bodyInfo;  // count, history head column, light delay on, history columns
    vec4 bodyShell;  // r range of all bodies, fraction of a step since the newest sample, step
};

const int kStar = 0, kRockyPlanet = 1, kGasGiant = 2, kPulsar = 3, kTidalStar = 4;

// Line segments of orbit trails and predicted paths, binned into a coarse grid.
struct PathSegment {
    vec4 a;      // start xyz, opacity
    vec4 b;      // end xyz, opacity
    vec4 color;  // rgb, dashed (> 0.5)
    vec4 dash;   // distance along the path at a and b
};
layout(std430, binding = 6) readonly buffer PathSegments {
    PathSegment segments[];
};
layout(std430, binding = 7) readonly buffer PathCells {
    uvec2 cellRange[];  // segments of cell c are pathIndex[cellRange[c].x .. + cellRange[c].y)
};
layout(std430, binding = 8) readonly buffer PathIndex {
    uint pathIndex[];
};

#include "noise.glsl"

const float PI = 3.14159265358979;

float spinA, spinA2;                // a, a^2
float photonE, photonL, photonQ;    // photon constants of motion (camera energy 1)

// --- Kerr geodesics -------------------------------------------------------------

// (P + vr) / Delta for the future-directed radial velocity vr (see kerr::inflowFactor):
// on the constraint surface it equals K / (P - vr), which is smooth where an
// ingoing ray crosses the horizon.
float inflowFactor(float P, float vr, float delta, float K) {
    if (P * vr < 0.0) return K / (P - vr);
    return (P + vr) / (abs(delta) > 1e-6 ? delta : (delta < 0.0 ? -1e-6 : 1e-6));
}

// Rates along the backward ray at (r, theta) with backward radial velocity vr.
void geodesicRhs(float r, float th, float vr, out float dPhi, out float accR, out float accTh, out float tDot) {
    float s = sin(th), c = cos(th);
    float s2 = max(s * s, 1e-8);
    float s3 = s2 * max(abs(s), 1e-4) * (s < 0.0 ? -1.0 : 1.0);
    float delta = r * r - 2.0 * r + spinA2;
    float rra = r * r + spinA2;
    float P = photonE * rra - spinA * photonL;
    float x = photonL - spinA * photonE;
    float K = photonQ + x * x;
    float F = inflowFactor(P, -vr, delta, K);
    accR = 2.0 * photonE * r * P - (r - 1.0) * K;
    accTh = -spinA2 * photonE * photonE * s * c + photonL * photonL * c / s3;
    dPhi = -(spinA * F - spinA * photonE + photonL / s2);
    tDot = rra * F - spinA * (spinA * photonE * s * s - photonL) + vr;
}

// Away from turning points, snap the velocities back onto vr^2 = R(r) and
// vtheta^2 = Theta(theta). This stops rounding errors made far from the hole,
// where R ~ r^4, from turning into fake turning points near it.
void enforceConstraints(float r, float th, inout float vr, inout float vth) {
    float delta = r * r - 2.0 * r + spinA2;
    float P = photonE * (r * r + spinA2) - spinA * photonL;
    float x = photonL - spinA * photonE;
    float K = photonQ + x * x;
    float R = P * P - delta * K;
    if (R > 0.05 * (P * P + abs(delta) * K)) vr = (vr < 0.0 ? -1.0 : 1.0) * sqrt(R);

    float s = sin(th), c = cos(th);
    float cot2 = c * c / max(s * s, 1e-8);
    float polar = spinA2 * photonE * photonE * c * c;
    float Theta = photonQ + polar - photonL * photonL * cot2;
    if (Theta > 0.05 * (abs(photonQ) + polar + photonL * photonL * cot2)) {
        vth = (vth < 0.0 ? -1.0 : 1.0) * sqrt(Theta);
    }
}

// Kerr-Schild Cartesian coordinates: x + i y = (r + i a) e^{i phi} sin(theta).
vec3 toCartesian(float r, float th, float ph) {
    float s = sin(th), cp = cos(ph), sp = sin(ph);
    return vec3((r * cp - spinA * sp) * s, (r * sp + spinA * cp) * s, r * cos(th));
}

vec3 cartesianVelocity(float r, float th, float ph, float vr, float vth, float dPh) {
    float s = sin(th), c = cos(th), sp = sin(ph), cp = cos(ph);
    float re = vr * s + r * c * vth - spinA * s * dPh;
    float im = spinA * c * vth + r * s * dPh;
    return vec3(re * cp - im * sp, re * sp + im * cp, vr * c - r * s * vth);
}

// Frequency ratio g = nu_camera / nu_emitter for gas or a body with 4-velocity u
// (Kerr-Schild components): 1 / (-p.u), the photon carrying unit camera energy.
// Covariant momentum: p_t = -E, p_phi = L, p_theta = -vth and p_r = F - E, from the
// backward-traced velocities vr and vth.
float frequencyRatio(vec4 u, float r, float vr, float vth) {
    float delta = r * r - 2.0 * r + spinA2;
    float P = photonE * (r * r + spinA2) - spinA * photonL;
    float x = photonL - spinA * photonE;
    float F = inflowFactor(P, -vr, delta, photonQ + x * x);
    float nuEmit = photonE * u.x - (F - photonE) * u.y + vth * u.z - photonL * u.w;
    return 1.0 / max(nuEmit, 1e-4);
}

// --- Color -----------------------------------------------------------------------

// Visible-band radiance of a blackbody at temperature T (luminance 1 at 6500 K).
vec3 blackbody(float T) {
    float n = float(textureSize(uBlackbody, 0).x);
    float u = (log2(max(T, 1.0)) * 0.30103 - uBlackbodyRange.x) / (uBlackbodyRange.y - uBlackbodyRange.x);
    u = clamp(u, 0.0, 1.0) * (n - 1.0) / n + 0.5 / n;
    return texture(uBlackbody, vec2(u, 0.5)).rgb;
}

vec3 blackbodyHue(float T) {
    vec3 c = blackbody(T);
    return c / max(max(c.r, max(c.g, c.b)), 1e-30);
}

float luminance(vec3 c) {
    return dot(c, vec3(0.2126, 0.7152, 0.0722));
}

// Applies the selected relativistic effects to a blackbody at rest-frame
// temperature T seen with frequency ratio g = nu_observed / nu_emitted.
// A blackbody stays a blackbody: it is simply seen at temperature g T, which
// includes the brightening (beaming) of approaching gas.
// The result is scaled by 1 / luminance(blackbody(Tref)), so that an emitter at
// the reference temperature has luminance 1 at rest whatever color it is.
vec3 shiftedBlackbody(float T, float g, float Tref) {
    float norm = 1.0 / max(luminance(blackbody(Tref)), 1e-30);
    if (uRelativity == 0) return blackbody(T * g) * norm;
    if (uRelativity == 1) return blackbody(T) * pow(g, 4.0) * norm;
    return blackbody(T) * norm;
}

// --- Accretion disk ---------------------------------------------------------------

float diskPattern(float r, float psi, float seed) {
    vec3 p = vec3(cos(psi) * 2.2, sin(psi) * 2.2, log(r) * 7.0 + seed);
    float broad = fbm(p);
    float fine = fbm(p * vec3(2.7, 2.7, 4.0) + 7.0);
    float rings = 0.5 + 0.5 * sin(log(r) * 38.0 + broad * 6.0);
    return clamp(broad * 0.6 + fine * 0.25 + rings * 0.2, 0.0, 1.0);
}

// Emission (rgb) and opacity (a) of the thin disk where a ray crosses it. delay is
// how long ago the light left that point: the turbulence is shown as it was then.
vec4 diskEmission(float r, float phi, float delay) {
    // Temperature profile of a thin disk with zero torque at the inner edge,
    // T^4 ~ r^-3 (1 - sqrt(r_in / r)), scaled so its peak equals uDiskTemp.
    float x = r / uDiskInner;
    float flux = max(pow(x, -3.0) * (1.0 - inversesqrt(x)), 0.0);
    float T = uDiskTemp * pow(flux / 0.0566524, 0.25);

    // Gas on circular equatorial geodesics.
    float omega = uDiskDirection / (pow(r, 1.5) + uDiskDirection * spinA);
    float gtt = -(1.0 - 2.0 / r);
    float gtp = -2.0 * spinA / r;
    float gpp = r * r + spinA2 + 2.0 * spinA2 / r;
    float norm = -(gtt + 2.0 * gtp * omega + gpp * omega * omega);
    if (norm <= 0.0) return vec4(0.0);  // no timelike circular orbit here
    float ut = inversesqrt(norm);
    float g = 1.0 / (ut * (photonE - photonL * omega));

    // Turbulence, sheared by differential rotation; two copies cross-fade so the
    // pattern is renewed before the shear winds it up completely.
    float phase1 = fract((uFlowTime - delay) / uFlowPeriod);
    float phase2 = fract(phase1 + 0.5);
    float n1 = diskPattern(r, phi - omega * phase1 * uFlowPeriod, 0.0);
    float n2 = diskPattern(r, phi - omega * phase2 * uFlowPeriod, 57.0);
    float n = mix(n1, n2, abs(2.0 * phase1 - 1.0));

    float density = mix(1.0, smoothstep(0.2, 0.75, n) * 1.5, uDiskTurbulence);
    float edges = (1.0 - smoothstep(uDiskOuter * 0.7, uDiskOuter, r)) *
                  smoothstep(uDiskInner, uDiskInner * 1.03, r);
    float alpha = clamp(uDiskOpacity * density * edges, 0.0, 1.0);

    vec3 emit = shiftedBlackbody(T, g, uDiskTemp) * uDiskBrightness * mix(1.0, 0.35 + 1.3 * n, uDiskTurbulence);
    return vec4(emit, alpha);
}

// --- Orbiting bodies ----------------------------------------------------------------

// The state of body i `delay` M of coordinate time ago, interpolated from the history.
struct BodySample {
    vec3 center;
    vec4 u;
    mat3 shape;  // symmetric: maps the body onto a sphere of its undisturbed radius
};

mat3 shapeMatrix(vec4 a, vec4 b) {
    return mat3(a.x, b.x, b.y, b.x, a.y, b.z, b.y, b.z, a.z);
}

// Which history samples describe body i `delay` M of coordinate time ago: columns c0
// and c1 and the fraction f between them (c1 < 0: between c0 and now). False if it
// has no history that far back, or had already gone by then.
struct HistoryRef {
    int row, c0, c1;
    float f;
};

bool historyRef(int i, float delay, out HistoryRef h) {
    h.row = kHistoryRows * int(bodies[i].motion.z + 0.5);
    h.c0 = -1;
    h.c1 = -1;
    h.f = 0.0;
    // A body that is gone is seen only in light that left it before it went.
    if (bodies[i].range.z >= 0.0 && delay <= bodies[i].range.z) return false;
    if (bodyInfo.z == 0 || delay <= 0.0) return true;  // as it is now
    float x = bodyShell.z - delay / bodyShell.w;  // in samples, relative to the newest one
    if (x >= 0.0) {
        h.c0 = bodyInfo.y;
        h.f = 1.0 - x / max(bodyShell.z, 1e-6);  // weight of the sample against "now"
        return true;
    }
    float i0 = floor(x);
    int back = int(-i0);
    if (back >= int(bodies[i].motion.w + 0.5)) return false;  // light from before its history
    h.c0 = (bodyInfo.y - back + bodyInfo.w) % bodyInfo.w;
    h.c1 = (h.c0 + 1) % bodyInfo.w;
    h.f = x - i0;
    return true;
}

// One history row of body i at that time (offset 0: position, 1: 4-velocity, 2 and 3:
// shape), interpolated; `now` is its current value.
vec4 historyRow(HistoryRef h, int offset, vec4 now) {
    if (h.c0 < 0) return now;
    vec4 a = texelFetch(uHistory, ivec2(h.c0, h.row + offset), 0);
    if (h.c1 < 0) return mix(now, a, h.f);
    return mix(a, texelFetch(uHistory, ivec2(h.c1, h.row + offset), 0), h.f);
}

// Per ray, for each body: how much farther the ray can go before it could possibly
// reach the body. Each step uses up its chord length, plus the distance the body can
// have moved meanwhile (its top speed times the change in light-travel time); while
// some is left, the body needs no look at all.
float bodySlack[32];

// Nearest hit of the segment p0 -> p1 with any body, as a fraction of the segment
// (2.0 for no hit). rLo..rHi is the r range the segment spans, delay the light-travel
// time to it and delayStep its change along the segment. Each body is an ellipsoid (a
// sphere unless tides deform it): mapped through its shape matrix, the test is
// against a sphere. Only bodies whose bounding sphere the segment passes through are
// looked at closely.
float hitBodies(vec3 p0, vec3 p1, float rLo, float rHi, float delay, float delayStep, out int which,
                out BodySample hit) {
    vec3 d = p1 - p0;
    float dd = max(dot(d, d), 1e-12);
    float chord = sqrt(dd);
    float best = 2.0;
    which = -1;
    for (int i = 0; i < bodyInfo.x; ++i) {
        bodySlack[i] -= chord + bodies[i].shapeB.w * abs(delayStep);
        if (bodySlack[i] > 0.0 || i == uSkipBody) continue;
        if (rHi < bodies[i].range.x || rLo > bodies[i].range.y) continue;
        HistoryRef h;
        if (!historyRef(i, delay, h)) continue;
        vec3 center = historyRow(h, 0, vec4(bodies[i].positionRadius.xyz, 0.0)).xyz;
        float bound = bodies[i].range.w;  // radius times the longest stretch
        bodySlack[i] = length(p1 - center) - bound;
        vec3 m0 = p0 - center;
        float along = clamp(-dot(m0, d) / dd, 0.0, 1.0);
        vec3 nearest = m0 + along * d;
        if (dot(nearest, nearest) > bound * bound) continue;

        mat3 shape = mat3(1.0);
        if (int(bodies[i].params.z + 0.5) == kTidalStar) {
            shape = shapeMatrix(historyRow(h, 2, bodies[i].shapeA), historyRow(h, 3, bodies[i].shapeB));
        }
        float radius = bodies[i].positionRadius.w;
        vec3 m = shape * m0;
        vec3 di = shape * d;
        float ddi = max(dot(di, di), 1e-12);
        float bq = dot(m, di);
        float cq = dot(m, m) - radius * radius;
        if (cq > 0.0 && bq > 0.0) continue;
        float disc = bq * bq - ddi * cq;
        if (disc < 0.0) continue;
        float t = cq < 0.0 ? 0.0 : (-bq - sqrt(disc)) / ddi;
        if (t <= 1.0 && t < best) {
            best = t;
            which = i;
            hit.center = center;
            hit.shape = shape;
            hit.u = historyRow(h, 1, bodies[i].fourVelocity);
        }
    }
    return best;
}

vec3 rockySurface(vec3 n, float seed) {
    float detail = fbm(n * 9.0 + seed * 1.7);
    float height = fbm(n * 2.5 + seed) - 0.5 + 0.25 * (detail - 0.5);
    vec3 ocean = mix(vec3(0.02, 0.06, 0.18), vec3(0.05, 0.18, 0.38), smoothstep(-0.2, 0.0, height));
    vec3 land = mix(vec3(0.13, 0.3, 0.09), vec3(0.45, 0.36, 0.22), smoothstep(0.03, 0.18, height));
    vec3 col = height < 0.0 ? ocean : land;
    col = mix(col, vec3(0.88, 0.92, 0.96), smoothstep(0.74, 0.8, abs(n.z) + 0.1 * detail));  // ice caps
    float clouds = smoothstep(0.55, 0.75, fbm(n * vec3(4.0, 4.0, 6.0) + seed * 3.1 + 11.0));
    return mix(col, vec3(0.92), clouds * 0.8);
}

vec3 gasGiantSurface(vec3 n, float seed) {
    float swirl = fbm(n * vec3(3.0, 3.0, 12.0) + seed);
    float bands = sin(n.z * 17.0 + swirl * 3.5);
    vec3 col = mix(vec3(0.5, 0.32, 0.2), vec3(0.88, 0.78, 0.62), 0.5 + 0.5 * bands);
    col *= 0.85 + 0.3 * fbm(n * 14.0 + seed * 2.0);
    vec3 stormDir = normalize(vec3(cos(seed), sin(seed), -0.35));
    float storm = smoothstep(0.93, 0.97, dot(normalize(n * vec3(1.0, 1.0, 2.2)), stormDir));
    return mix(col, vec3(0.78, 0.32, 0.2), storm);
}

// Light leaving body i toward the camera. rayDir is the backward ray's direction,
// so -rayDir is the direction the light travels; delay is how long ago it left.
vec3 bodyRadiance(int i, vec3 hit, BodySample b, vec3 rayDir, float r, float vr, float vth, float delay) {
    // The ellipsoid's normal is the gradient of |S x|^2, S^2 x; the surface pattern
    // is painted on the undeformed sphere, S x.
    vec3 offset = hit - b.center;
    vec3 sphere = b.shape * offset;
    vec3 normal = normalize(b.shape * sphere);
    vec3 center = b.center;
    vec3 toCamera = -rayDir;
    float g = frequencyRatio(b.u, r, vr, vth);
    int kind = int(bodies[i].params.z + 0.5);
    float seed = bodies[i].params.w;
    float brightness = bodies[i].params.y;
    // The body's rotation when the light left, and the surface normal in its own frame.
    float angle = bodies[i].motion.x - bodies[i].motion.y * delay;
    float ca = cos(angle), sa = sin(angle);
    vec3 n0 = normalize(sphere);
    vec3 local = vec3(ca * n0.x + sa * n0.y, -sa * n0.x + ca * n0.y, n0.z);

    if (kind == kRockyPlanet || kind == kGasGiant) {
        // Lit by the hot inner disk, treated as a light at the center falling off as 1/d^2.
        vec3 albedo = kind == kRockyPlanet ? rockySurface(local, seed) : gasGiantSurface(local, seed);
        float irradiance = uDiskEnabled == 1 ? uDiskBrightness * 36.0 / max(dot(center, center), 4.0) : 0.0;
        float lambert = max(dot(normal, normalize(-center)), 0.0);
        return albedo * shiftedBlackbody(uDiskTemp, g, uDiskTemp) * (irradiance * lambert + 0.01) * brightness;
    }

    float T = bodies[i].params.x;
    float limb = 0.4 + 0.6 * clamp(dot(normal, toCamera), 0.0, 1.0);  // linear limb darkening
    vec3 color = shiftedBlackbody(T, g, T) * brightness * limb;
    if (kind == kPulsar) {
        // Two beams along a magnetic axis tilted 50 degrees from the rotation axis.
        // We see a flash when one points along the light's direction of travel.
        vec3 axis = vec3(sin(0.87) * cos(angle), sin(0.87) * sin(angle), cos(0.87));
        float beam = exp(-(1.0 - abs(dot(toCamera, axis))) / 0.012);
        return color * (1.5 + 60.0 * beam);
    }
    return color * (0.8 + 0.4 * fbm(local * 7.0 + seed));  // granulation
}

// --- Debris gas ----------------------------------------------------------------------

vec3 gridCoord(float r, float th, float ph) {
    return vec3((log(r) - uGridLogR.x) / (uGridLogR.y - uGridLogR.x), th / PI, ph / (2.0 * PI));
}

bool debrisNear(vec4 p) {
    ivec3 dims = textureSize(uDebrisOccupancy, 0);
    vec3 g = gridCoord(p.x, p.y, p.z);
    ivec3 cell = ivec3(floor(vec3(g.x, g.y, fract(g.z)) * vec3(dims)));
    return texelFetch(uDebrisOccupancy, clamp(cell, ivec3(0), dims - 1), 0).x != 0u;
}

// The averaged 4-velocity of the gas in a voxel is not quite a unit vector (it mixes
// parcels, and is stored at half precision): rescale it to u.u = -1 with the
// Kerr-Schild metric. Returns false if it is not timelike.
bool normalizeVelocity(inout vec4 u, float r, float th) {
    float s = sin(th), c = cos(th);
    float s2 = s * s;
    float sigma = r * r + spinA2 * c * c;
    float f = 2.0 * r / sigma;
    float rra = r * r + spinA2;
    float bigA = rra * rra - spinA2 * (r * r - 2.0 * r + spinA2) * s2;
    float uu = -(1.0 - f) * u.x * u.x + 2.0 * f * u.x * u.y - 2.0 * f * spinA * s2 * u.x * u.w +
               (1.0 + f) * u.y * u.y - 2.0 * spinA * s2 * (1.0 + f) * u.y * u.w + sigma * u.z * u.z +
               bigA * s2 / sigma * u.w * u.w;
    if (uu > -1e-4) return false;
    u *= inversesqrt(-uu);
    return true;
}

// Light from hot, thin gas at temperature T seen with frequency ratio g: a blackbody
// color, with brightness following the heat it radiates (~T), boosted by g^3 for
// optically thin emission.
vec3 gasGlow(float T, float g) {
    vec3 hue = blackbody(uRelativity == 0 ? T * g : T);
    hue /= max(luminance(hue), 1e-30);
    float beaming = uRelativity == 2 ? 1.0 : g * g * g;
    return hue * (T / 6000.0) * beaming;
}

// Emission and absorption by debris along one ray segment (two endpoints with their
// coordinates and backward-traced velocities). With light delay on, the grid holds
// each parcel where it was when the light now reaching the camera left it.
void traceDebris(vec4 a, vec2 va, vec4 b, vec2 vb, float segLength, inout vec3 color, inout float transmittance) {
    float rMax = exp(uGridLogR.y);
    if (min(a.x, b.x) > rMax) return;
    // Skip empty space using the coarse occupancy map.
    if (!debrisNear(a) && !debrisNear(b) && !debrisNear(0.5 * (a + b))) return;

    float voxel = 0.025 * min(a.x, b.x);  // grid spacing at this radius
    int n = clamp(int(ceil(segLength / (0.5 * voxel))), 1, 16);
    float dl = segLength / float(n);
    for (int k = 0; k < n; ++k) {
        float s = (float(k) + 0.5) / float(n);
        vec4 p = mix(a, b, s);  // r, theta, phi
        if (p.x > rMax) continue;
        vec3 uvw = gridCoord(p.x, p.y, p.z);
        vec4 d0 = textureLod(uDebrisDensity, uvw, 0.0);
        float rho = d0.x;
        if (rho <= 1e-5) continue;
        vec2 d1 = textureLod(uDebrisVelocity, uvw, 0.0).xy;
        vec4 u = vec4(d0.z, d0.w, d1.x, d1.y) / rho;
        if (!normalizeVelocity(u, p.x, p.y)) continue;
        float T = d0.y / rho * 1e4;
        vec2 v = mix(va, vb, s);
        // Beaming as strong as g = 4 needs gas moving at 0.9c straight at us; anything
        // more comes from mixing parcels moving in different directions.
        float g = min(frequencyRatio(u, p.x, v.x, v.y), 4.0);
        // Emission (glow) and absorption (thickness) are set separately.
        color += transmittance * gasGlow(T, g) * (uDebrisBrightness * rho * dl);
        transmittance *= exp(-rho * uDebrisOpacity * dl);
        if (transmittance < 0.002) return;
    }
}

// --- Orbit paths --------------------------------------------------------------------

// The path cell holding p (grid coordinates), and whether it holds any segments (a
// cheap texture lookup that spares most steps the reads of the segment lists).
bool pathCell(vec4 p, out ivec3 c) {
    vec3 g = gridCoord(p.x, p.y, p.z);
    c = ivec3(-1000);
    if (g.x < 0.0 || g.x >= 1.0) return false;
    c = min(ivec3(floor(vec3(g.x, clamp(g.y, 0.0, 0.9999), fract(g.z)) * vec3(uPathGridDims))), uPathGridDims - 1);
    return texelFetch(uPathOccupancy, c, 0).x != 0u;
}

int pathCellIndex(ivec3 c) {
    return (c.z * uPathGridDims.y + c.y) * uPathGridDims.x + c.x;
}

// Glow of the path lines of one cell passing close to the ray chord x0 -> x1.
// pathLength is how far the light has travelled from the camera to x0, so lines keep
// a constant width on screen. Images seen through strongly bent light (around the
// hole, near the photon ring) are squeezed far thinner than a pixel in reality; they
// fade with the bending, measured as straight-line distance over path length.
// Segments already drawn for this chord (from another cell) are in `seen`.
void tracePathCell(int cell, vec3 x0, vec3 x1, float pathLength, vec3 cameraX, inout vec3 color, float transmittance,
                   inout uint seen[4], inout int seenCount) {
    vec3 d1 = x1 - x0;
    float a11 = dot(d1, d1);
    float chord = sqrt(a11);
    uvec2 range = cellRange[cell];
    for (uint k = range.x; k < range.x + range.y; ++k) {
        uint index = pathIndex[k];
        vec3 sa = segments[index].a.xyz;
        vec3 d2 = segments[index].b.xyz - sa;
        // Closest points between the chord (parameter u) and the segment (v).
        vec3 w0 = x0 - sa;
        float a22 = dot(d2, d2), a12 = dot(d1, d2);
        float b1 = dot(d1, w0), b2 = dot(d2, w0);
        float den = a11 * a22 - a12 * a12;
        float u = den > 1e-12 ? clamp((a12 * b2 - a22 * b1) / den, 0.0, 1.0) : 0.0;
        float v = a22 > 1e-12 ? clamp((a12 * u + b2) / a22, 0.0, 1.0) : 0.0;
        u = a11 > 1e-12 ? clamp((a12 * v - b1) / a11, 0.0, 1.0) : 0.0;
        // Each crossing counts once: skip shared end points of chords and segments.
        if (u >= 1.0 || v >= 1.0) continue;
        float width = uPathWidth * uPixelAngle * (pathLength + u * chord);
        vec3 gap = (x0 + u * d1) - (sa + v * d2);
        float dist2 = dot(gap, gap);
        if (dist2 > 9.0 * width * width) continue;
        bool duplicate = false;
        for (int j = 0; j < seenCount; ++j) duplicate = duplicate || seen[j] == index;
        if (duplicate) continue;
        if (seenCount < 4) seen[seenCount++] = index;
        PathSegment s = segments[index];
        float along = mix(s.dash.x, s.dash.y, v);
        float dash = s.color.w > 0.5 ? step(0.45, fract(along * 0.8)) : 1.0;
        float opacity = mix(s.a.w, s.b.w, v);
        float straightness = clamp(length(x0 + u * d1 - cameraX) / max(pathLength + u * chord, 1e-3), 0.0, 1.0);
        float bentFade = straightness * straightness * straightness;
        color += transmittance * s.color.rgb * (opacity * dash * bentFade * exp(-dist2 / (width * width)));
    }
}

// Path lines near the chord x0 -> x1 (coordinates p0, p1). c0 is p0's cell, carried
// over from the previous step; c1 is set to p1's. Cells are checked only if they
// hold segments; the midpoint's cell only when the chord crosses cells diagonally.
void tracePaths(vec3 x0, vec3 x1, vec4 p0, vec4 p1, float pathLength, vec3 cameraX, ivec3 c0, bool occupied0,
                out ivec3 c1, out bool occupied1, inout vec3 color, float transmittance) {
    occupied1 = pathCell(p1, c1);
    uint seen[4] = uint[4](0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu);
    int seenCount = 0;
    if (occupied0) tracePathCell(pathCellIndex(c0), x0, x1, pathLength, cameraX, color, transmittance, seen, seenCount);
    if (occupied1 && c1 != c0) {
        tracePathCell(pathCellIndex(c1), x0, x1, pathLength, cameraX, color, transmittance, seen, seenCount);
    }
    ivec3 jump = abs(c1 - c0);
    jump.z = min(jump.z, uPathGridDims.z - jump.z);
    if (c0.x > -1000 && c1.x > -1000 && (jump.x + jump.y + jump.z > 1)) {
        ivec3 cm;
        if (pathCell(0.5 * (p0 + p1), cm) && cm != c0 && cm != c1) {
            tracePathCell(pathCellIndex(cm), x0, x1, pathLength, cameraX, color, transmittance, seen, seenCount);
        }
    }
}

// --- Relativistic jets -------------------------------------------------------------------

// Twin jets of plasma along the spin axis, moving at uJetSpeed. The approaching jet
// is boosted by Doppler beaming; the receding one fades. Emission ~ g^3, as for
// optically thin synchrotron light with a flat spectrum.
// delay0 and delay1 are the light-travel times to the chord's ends: the knots are
// drawn where they were when the light left them, which makes a jet coming toward
// us appear to move faster than light.
void traceJets(vec3 x0, vec3 x1, vec4 p0, vec4 p1, vec2 v0, vec2 v1, float delay0, float delay1, float horizon,
               inout vec3 color, float transmittance) {
    const float zMax = 90.0;
    float z0 = horizon + 1.0;
    float zm = max(abs(x0.z), abs(x1.z));
    if (min(abs(x0.z), abs(x1.z)) > zMax || zm < z0) return;
    float wMax = 0.4 + zm * uJetOpening;
    // Skip chords that stay far from the axis.
    vec2 a = x0.xy, d = x1.xy - x0.xy;
    float t = clamp(-dot(a, d) / max(dot(d, d), 1e-12), 0.0, 1.0);
    if (length(a + t * d) > 3.0 * wMax) return;

    float chord = length(x1 - x0);
    int n = clamp(int(ceil(chord / (0.4 * (0.4 + 0.5 * zm * uJetOpening)))), 1, 12);
    float dl = chord / float(n);
    vec3 jetColor = vec3(0.55, 0.7, 1.0);
    for (int k = 0; k < n; ++k) {
        float s = (float(k) + 0.5) / float(n);
        vec3 x = mix(x0, x1, s);
        float z = abs(x.z);
        if (z < z0 || z > zMax) continue;
        float w = 0.4 + z * uJetOpening;
        float rhoCyl = length(x.xy);
        float profile = exp(-(rhoCyl * rhoCyl) / (w * w)) * pow(z0 / z, 1.2) * smoothstep(z0, z0 + 2.0, z) *
                        (1.0 - smoothstep(60.0, zMax, z));
        if (profile < 1e-4) continue;
        // Knots of brighter plasma moving out along the jet.
        float emitted = uTime - mix(delay0, delay1, s);
        float knots = 0.55 + 0.45 * sin(z * 0.45 - emitted * uJetSpeed * 0.45);
        profile *= knots * knots;

        // Plasma 4-velocity: speed beta along the axis, away from the hole, in the ZAMO frame.
        vec4 p = mix(p0, p1, s);
        float th = p.y;
        float sth = max(sin(th), 1e-4), cth = cos(th);
        float side = x.z >= 0.0 ? 1.0 : -1.0;
        vec2 vLocal = side * uJetSpeed * vec2(cth, -sth);  // (r, theta) components of +-z
        float r = p.x;
        float sigma = r * r + spinA2 * cth * cth;
        float delta = r * r - 2.0 * r + spinA2;
        float rra = r * r + spinA2;
        float bigA = rra * rra - spinA2 * delta * sth * sth;
        float gamma = inversesqrt(1.0 - uJetSpeed * uJetSpeed);
        float ut = gamma * sqrt(bigA / (sigma * delta));
        float ur = gamma * vLocal.x * sqrt(delta / sigma);
        // Boyer-Lindquist components, shifted to Kerr-Schild ones.
        vec4 u = vec4(ut + (rra / delta - 1.0) * ur, ur, gamma * vLocal.y / sqrt(sigma),
                      ut * 2.0 * spinA * r / bigA + spinA / delta * ur);
        vec2 v = mix(v0, v1, s);
        float g = frequencyRatio(u, r, v.x, v.y);
        float boost = uRelativity == 2 ? 1.0 : g * g * g;
        color += transmittance * jetColor * (uJetPower * 0.6 * profile * boost * dl);
    }
}

// --- Sky ---------------------------------------------------------------------------

vec3 starLayer(vec3 d, float cells, float density, float seed, float magnitude, float gShift) {
    vec3 p = d * cells;
    vec3 base = floor(p);
    // Keep each star well inside the 3x3x3 cells searched, or it is clipped square.
    float sigma = min(max(uPixelAngle * 0.7, 2e-4), 0.3 / cells);
    vec3 col = vec3(0.0);
    for (int z = -1; z <= 1; ++z)
        for (int y = -1; y <= 1; ++y)
            for (int x = -1; x <= 1; ++x) {
                vec3 cell = base + vec3(x, y, z);
                vec4 h = hash43(cell + seed);
                if (h.w > density) continue;
                vec3 star = normalize(cell + 0.5 + (h.xyz - 0.5) * 0.8);
                vec3 diff = d - star;
                float falloff = exp(-dot(diff, diff) / (2.0 * sigma * sigma));
                if (falloff < 1e-4) continue;
                vec4 k = hash43(cell + seed + 71.3);
                float bright = magnitude * (0.015 + 5.0 * pow(k.x, 14.0));
                float temp = mix(2800.0, 15000.0, k.y * k.y);
                col += blackbodyHue(temp * gShift) * bright * falloff;
            }
    return col;
}

vec3 galaxy(vec3 d) {
    const vec3 pole = normalize(vec3(-0.35, 0.55, 0.76));
    const vec3 core = normalize(vec3(0.85, 0.42, 0.05));
    float lat = asin(clamp(dot(d, pole), -1.0, 1.0));
    float band = exp(-lat * lat / (2.0 * 0.16 * 0.16));
    float angle = acos(clamp(dot(d, core), -1.0, 1.0));
    float bulge = exp(-angle * angle / (2.0 * 0.4 * 0.4));
    float clouds = fbm(d * 4.0 + 3.1);
    float dust = smoothstep(0.45, 0.75, fbm(d * 9.0 + 11.7)) * exp(-lat * lat / (2.0 * 0.05 * 0.05));
    vec3 tint = mix(vec3(0.45, 0.55, 0.9), vec3(1.0, 0.8, 0.6), clamp(bulge * 1.5 + clouds * 0.4, 0.0, 1.0));
    float intensity = band * (0.3 + clouds) * (1.0 + 2.5 * bulge) * (1.0 - 0.8 * dust);
    return tint * intensity * 0.05 * uGalaxyBrightness;
}

vec3 gridSky(vec3 d) {
    float lon = atan(d.y, d.x);
    float lat = asin(clamp(d.z, -1.0, 1.0));
    vec2 uv = vec2(lon, lat) / PI * 12.0;  // 15-degree squares
    vec2 cell = floor(uv);
    int quadrant = int(floor((lon + PI) / (0.5 * PI))) & 3;
    vec3 base = quadrant == 0 ? vec3(0.95, 0.3, 0.25)
              : quadrant == 1 ? vec3(0.25, 0.75, 0.35)
              : quadrant == 2 ? vec3(0.3, 0.45, 0.95)
                              : vec3(0.95, 0.8, 0.25);
    if (lat < 0.0) base *= 0.55;
    vec3 col = base * mix(0.45, 1.0, mod(cell.x + cell.y, 2.0));
    vec2 edge = 0.5 - abs(fract(uv) - 0.5);
    vec2 lineWidth = max(fwidth(uv), vec2(1e-4)) * 1.2;
    vec2 line = smoothstep(vec2(0.0), lineWidth, edge);
    return col * mix(0.15, 1.0, min(line.x, line.y)) * 0.6;
}

vec3 imageSky(vec3 d) {
    float lon = atan(d.y, d.x);
    float lat = asin(clamp(d.z, -1.0, 1.0));
    vec2 uv = vec2(lon / (2.0 * PI) + 0.5, lat / PI + 0.5);
    return textureLod(uSkyTex, uv, uSkyLod).rgb * uSkyImageBrightness;
}

// gSky: the frequency ratio for light from far away, 1 / E.
vec3 skyColor(vec3 d, float gSky) {
    if (uSkyMode == 1) return gridSky(d);
    if (uSkyMode == 2) return imageSky(d);
    if (uSkyMode == 3) return vec3(0.0);
    float shift = uRelativity == 0 ? gSky : 1.0;  // starlight is blueshifted falling toward the hole
    vec3 col = galaxy(d);
    col += starLayer(d, 40.0, 0.25, 0.0, 1.0, shift);
    col += starLayer(d, 100.0, 0.06, 17.0, 0.6, shift);
    col += starLayer(d, 220.0, 0.025, 43.0, 0.35, shift);
    return col * uStarBrightness;
}

vec3 heatMap(float x) {
    return clamp(vec3(1.5 * x, 1.5 * x * (1.0 - x) * 2.0, 1.0 - 1.5 * x), 0.0, 1.0);
}

// --- Main ----------------------------------------------------------------------------

void main() {
    spinA = uSpin;
    spinA2 = spinA * spinA;

    // Direction of this pixel in the camera's rest frame, components (r, theta, phi).
    vec2 ndc = (gl_FragCoord.xy - 0.5 * uResolution) / (0.5 * uResolution.y);
    vec3 n = normalize(uCamBasis * vec3(ndc * uTanHalfFov, 1.0));

    // The photon arriving from n, with unit energy in the camera frame (see
    // kerr::photonFromFrame). A moving camera's frame is boosted, so the view crowds
    // toward the direction of motion and is blueshifted there: aberration.
    float r = uCamPos.x, th = uCamPos.y, ph = uCamPos.z;
    vec4 q = uCamPhoton * vec4(1.0, -n);
    photonE = -q.x;
    photonL = q.y;
    float pTheta = q.z;
    float s = sin(th), c = cos(th);
    float sigma = r * r + spinA2 * c * c;
    photonQ = pTheta * pTheta + c * c * (photonL * photonL / max(s * s, 1e-8) - spinA2 * photonE * photonE);
    // Backward-traced velocities: the photon moves toward the camera.
    float vr = -sigma * q.w;
    float vth = -pTheta;
    if (abs(photonL) < 3e-4 * sqrt(max(photonQ + spinA2 * photonE * photonE, 0.0))) {
        // Would pass within ~3e-4 rad of the spin axis, where phi swings through pi
        // too fast to resolve; send it exactly through the axis (kerr::kAxisSnap).
        photonL = 0.0;
        photonQ = pTheta * pTheta - c * c * spinA2 * photonE * photonE;
    }
    float gSky = photonE > 1e-6 ? 1.0 / photonE : 1.0;

    float horizon = 1.0 + sqrt(max(1.0 - spinA2, 0.0));
    bool lightDelay = bodyInfo.z != 0;
    vec3 X = toCartesian(r, th, ph);
    vec3 cameraX = X;
    vec3 color = vec3(0.0);
    float transmittance = 1.0;
    float delay = 0.0;       // coordinate time since the light was at the current point
    float pathLength = 0.0;  // distance the light travelled from the current point to the camera
    int fate = 0;       // 0 = undecided, 1 = escaped to the sky, 2 = captured or blocked
    vec3 escapeDir = vec3(0.0, 0.0, 1.0);
    int steps = 0;

    // Inside the horizon every ray, followed back in time, rises to the horizon. If it
    // gets there with P(r+) = 2 r+ E - a L < 0 it entered from the far side of the
    // eternal Kerr spacetime, from which no light is drawn: dark.
    if (r < horizon && 2.0 * horizon * photonE - spinA * photonL < 0.0) fate = 2;

    for (int k = 0; k < 32; ++k) bodySlack[k] = -1.0;
    ivec3 pathCellNow = ivec3(-1000);
    bool pathOccupiedNow = uPathsEnabled == 1 && pathCell(vec4(r, th, ph, 0.0), pathCellNow);

    for (int i = 0; i < uMaxSteps && fate == 0; ++i) {
        steps = i + 1;
        float dPh1, aR1, aTh1, tD1;
        geodesicRhs(r, th, vr, dPh1, aR1, aTh1, tD1);

        // Adaptive step (kerr::photonStepSize is the tested CPU twin): bound the
        // relative change of r, of sin(theta) and of phi. Near the spin axis phi
        // swings through ~pi within a tiny distance, so steps must shrink as the
        // ray approaches the axis, before the swing starts. Far from the hole
        // spacetime is nearly flat, so the allowed change grows with r.
        float sn = abs(sin(th));
        float tol = uStepSize * clamp(r / 8.0, 1.0, 4.0);
        float h = tol / (abs(vr) / r + abs(vth) / max(sn, 1e-4) + abs(dPh1) + 1e-6);
        h = min(h, sqrt(tol * r / (abs(aR1) + 1e-6)));
        h = min(h, sqrt(tol * max(sn, 1e-4) / (abs(aTh1) + 1e-6)));

        // Classical RK4 on (r, theta, phi, vr, vtheta), and on t for the delay.
        float dPh2, aR2, aTh2, tD2, dPh3, aR3, aTh3, tD3, dPh4, aR4, aTh4, tD4;
        float vr2 = vr + 0.5 * h * aR1, vth2 = vth + 0.5 * h * aTh1;
        geodesicRhs(r + 0.5 * h * vr, th + 0.5 * h * vth, vr2, dPh2, aR2, aTh2, tD2);
        float vr3 = vr + 0.5 * h * aR2, vth3 = vth + 0.5 * h * aTh2;
        geodesicRhs(r + 0.5 * h * vr2, th + 0.5 * h * vth2, vr3, dPh3, aR3, aTh3, tD3);
        float vr4 = vr + h * aR3, vth4 = vth + h * aTh3;
        geodesicRhs(r + h * vr3, th + h * vth3, vr4, dPh4, aR4, aTh4, tD4);

        float r1 = r + h / 6.0 * (vr + 2.0 * vr2 + 2.0 * vr3 + vr4);
        float th1 = th + h / 6.0 * (vth + 2.0 * vth2 + 2.0 * vth3 + vth4);
        float ph1 = ph + h / 6.0 * (dPh1 + 2.0 * dPh2 + 2.0 * dPh3 + dPh4);
        float vrNext = vr + h / 6.0 * (aR1 + 2.0 * aR2 + 2.0 * aR3 + aR4);
        float vthNext = vth + h / 6.0 * (aTh1 + 2.0 * aTh2 + 2.0 * aTh3 + aTh4);
        float delayNext = delay + h / 6.0 * (tD1 + 2.0 * tD2 + 2.0 * tD3 + tD4);

        // Crossing a pole: same point, described with theta back in [0, pi].
        bool reflected = false;
        if (th1 < 0.0) {
            th1 = -th1;
            ph1 += PI;
            vthNext = -vthNext;
            reflected = true;
        } else if (th1 > PI) {
            th1 = 2.0 * PI - th1;
            ph1 += PI;
            vthNext = -vthNext;
            reflected = true;
        }
        enforceConstraints(r1, th1, vrNext, vthNext);
        vec3 X1 = toCartesian(r1, th1, ph1);

        // Glowing matter along the whole segment: debris gas, jets and path lines
        // (not across a pole, where the coordinates jump).
        float chordLength = length(X1 - X);
        vec4 p0 = vec4(r, th, ph, 0.0), p1 = vec4(r1, th1, ph1, 0.0);
        if (uDebrisEnabled == 1 && !reflected) {
            traceDebris(p0, vec2(vr, vth), p1, vec2(vrNext, vthNext), chordLength, color, transmittance);
        }
        if (uJetsEnabled == 1 && !reflected) {
            traceJets(X, X1, p0, p1, vec2(vr, vth), vec2(vrNext, vthNext), lightDelay ? delay : 0.0,
                      lightDelay ? delayNext : 0.0, horizon, color, transmittance);
        }
        if (uPathsEnabled == 1) {
            ivec3 nextCell;
            bool nextOccupied;
            if (!reflected) {
                tracePaths(X, X1, p0, p1, pathLength, cameraX, pathCellNow, pathOccupiedNow, nextCell, nextOccupied,
                           color, transmittance);
            } else {
                nextOccupied = pathCell(p1, nextCell);
            }
            pathCellNow = nextCell;
            pathOccupiedNow = nextOccupied;
        }
        pathLength += chordLength;

        // Events along this segment, handled in order: disk crossing, body hit.
        float tDisk = 2.0;
        // Strictly across: a camera sitting exactly in the plane is not on the disk.
        if (uDiskEnabled == 1 && (th - 0.5 * PI) * (th1 - 0.5 * PI) < 0.0) {
            tDisk = (th - 0.5 * PI) / (th - th1);
        }
        int body = -1;
        float tBody = 2.0;
        BodySample bodyHit;
        float rLo = min(r, r1), rHi = max(r, r1);
        float delayMid = lightDelay ? 0.5 * (delay + delayNext) : 0.0;
        if (bodyInfo.x > 0 && rHi > bodyShell.x && rLo < bodyShell.y) {
            tBody = hitBodies(X, X1, rLo, rHi, delayMid, lightDelay ? delayNext - delay : 0.0, body, bodyHit);
        } else {
            for (int k = 0; k < bodyInfo.x; ++k) bodySlack[k] -= chordLength;
        }
        if (tDisk <= 1.0 && tDisk < tBody) {
            float rc = mix(r, r1, tDisk);
            if (rc >= uDiskInner && rc <= uDiskOuter) {
                float diskDelay = lightDelay ? mix(delay, delayNext, tDisk) : 0.0;
                vec4 e = diskEmission(rc, mix(ph, ph1, tDisk), diskDelay);
                color += transmittance * e.a * e.rgb;
                transmittance *= 1.0 - e.a;
            }
        }
        if (tBody <= 1.0) {
            vec3 hit = mix(X, X1, tBody);
            color += transmittance * bodyRadiance(body, hit, bodyHit, normalize(X1 - X),
                                                  mix(r, r1, tBody), mix(vr, vrNext, tBody),
                                                  mix(vth, vthNext, tBody), delayMid);
            transmittance = 0.0;
        }

        r = r1;
        th = th1;
        ph = ph1;
        vr = vrNext;
        vth = vthNext;
        delay = delayNext;
        X = X1;

        // Sinking back toward the horizon from outside: the light came out of the past
        // horizon, which shows as the black hole's shadow.
        if (transmittance < 0.002 || (r < horizon + 0.01 && vr < 0.0)) {
            fate = 2;
            break;
        }
        if (r > uEscapeRadius && vr > 0.0) {
            float dPh, unusedR, unusedTh, unusedT;
            geodesicRhs(r, th, vr, dPh, unusedR, unusedTh, unusedT);
            escapeDir = normalize(cartesianVelocity(r, th, ph, vr, vth, dPh));
            fate = 1;
            break;
        }
    }
    if (fate == 0 && r > 10.0) {
        // Out of steps far from the hole: close enough to escaping.
        escapeDir = normalize(X);
        fate = 1;
    }

    // Evaluated outside the loop so screen-space derivatives (grid lines) are valid.
    vec3 sky = skyColor(escapeDir, gSky);
    if (fate == 1) color += transmittance * sky;

    if (uDebugView == 1) color = heatMap(float(steps) / float(uMaxSteps));
    if (uDebugView == 2) color = fate == 1 ? escapeDir * 0.5 + 0.5 : vec3(0.0);
    if (any(isnan(color)) || any(isinf(color))) color = vec3(0.0);
    outColor = vec4(min(color, vec3(30000.0)), 1.0);
}
