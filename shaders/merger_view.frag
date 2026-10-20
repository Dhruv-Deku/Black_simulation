#version 460 core

// The neutron star merger, cinematic view (src/merger_view.*).
//
// Each pixel's ray is marched from the camera through the scene, bent by the stars' gravity:
// outside a spherical mass M a light ray obeys u'' + u = 3 M u^2 (u = 1 / r), which in
// Cartesian form is the acceleration
//     d^2 x / ds^2 = -3 M h^2 x / r^5,   h = |x  x  dx/ds|,
// exact for one Schwarzschild mass (the photon sphere at 3M and the shadow of radius
// sqrt(27) M come out of it). Two stars are treated as two such lenses added together
// (an approximation: the spacetime of a merging binary has no closed form), the black hole
// afterwards as one that swallows what crosses 2M.
//
// Along the ray:
//   the stars and the remnant  the simulation's matter (a 3D grid): it absorbs and glows by
//                              its density and temperature, Doppler-boosted by its motion,
//                              and textured by where each parcel started in its star, so
//                              the patterns stretch with the matter as tides pull it apart;
//   gravitational waves        drawn (in false color) where the ray crosses the orbital
//                              plane: the quadrupole pattern Re[H(t - r) e^(-2 i phi)] / r,
//                              the two-armed spiral of waves leaving at the speed of light;
//   the debris cloud           what the collision throws out (src/merger.h, DebrisCloud):
//                              two trailing tidal tails in the orbital plane and the shocked
//                              gas puffing out all round, billowing, hot inside and cooling
//                              as it expands, lit by the remnant within;
//   the kilonova               three homologously expanding components (src/merger.h),
//                              each glowing at its photosphere temperature, its optical
//                              depth from its opacity, mass and size;
//   the jet and cocoon         a narrow jet drilling out along the spin axis, and after it
//                              breaks out the hot cocoon spilling sideways.
// Units: render units (src/merger_view.cpp chooses them to keep numbers near 1).

layout(location = 0) out vec4 outColor;

#include "noise.glsl"

uniform vec2 uResolution;
uniform vec2 uViewport;       // lower-left pixel of this view in the target
uniform vec3 uCamPos;
uniform mat3 uCamBasis;       // right, up, forward
uniform float uTanHalfFov;
uniform float uPixelAngle;
uniform float uTime;          // seconds of physical time since the start (for animated detail)
uniform int uFilm;            // film colors

// Lenses: position, mass (render units). uHole: index of the black hole (-1: none).
uniform vec4 uLens[2];
uniform int uLensCount;
uniform int uHole;

// The simulation's grid.
layout(binding = 1) uniform sampler3D uVolumeA;  // log10 rho (g/cm^3), log10 T, label xy
layout(binding = 2) uniform sampler3D uVolumeB;  // velocity (c), label z
uniform int uVolumeOn;
uniform vec3 uBoxMin;
uniform float uBoxSize;
uniform float uOpacityScale;   // extinction per render unit at 1e15 g/cm^3
uniform float uMatterGlow;
uniform float uLabelScale;     // label units (stellar radii) per render unit

// Gravitational waves: face-on complex strain H(t) (normalized), sampled in time.
layout(binding = 3) uniform sampler2D uWaves;
uniform int uWavesOn;
uniform float uWaveNow;        // the texture's coordinate for now
uniform float uWaveStep;       // texture coordinate per render unit of distance (light travel)
uniform float uWaveInner, uWaveOuter;
uniform float uWaveGlow;

// The debris cloud: edge radius, optical depth through its radius (as drawn), brightness,
// temperature (K, as shown); the remnant's radius (the cloud starts outside it), the
// stars' line of contact (angle in the orbital plane), how far the billows have churned,
// how brightly the remnant lights it from within.
uniform vec4 uCloud;
uniform vec4 uCloudMore;
uniform float uMergeFlash;    // the collision's flash
// Tileable noise for the clouds (merger_noise.comp): billows, ridges, two to warp with.
layout(binding = 5) uniform sampler3D uNoise;

// Kilonova: per component outer radius, display temperature, brightness, optical depth
// through its radius.
uniform vec4 uKilonova[3];
uniform int uKilonovaOn;
layout(binding = 4) uniform sampler2D uEjectaSky;  // where the dynamical ejecta went (from the simulation)

// Jet: head distance along the axis, cocoon radius, brightness, half-opening angle;
// breakout distance, afterglow (X-ray) brightness.
uniform vec4 uJet;
uniform vec2 uJetMore;
uniform vec3 uJetAxis;
uniform float uFlash;          // the gamma-ray flash

layout(binding = 0) uniform sampler2D uBlackbody;
uniform vec2 uBlackbodyRange;
uniform float uStarBrightness;
uniform float uGalaxyBrightness;

const float PI = 3.14159265358979;

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

// --- Sky (as in trace.frag) ------------------------------------------------------------

vec3 starLayer(vec3 d, float cells, float density, float seed, float magnitude) {
    vec3 p = d * cells;
    vec3 base = floor(p);
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
                col += blackbodyHue(temp) * bright * falloff;
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

vec3 sky(vec3 d) {
    vec3 col = galaxy(d);
    col += starLayer(d, 40.0, 0.25, 0.0, 1.0);
    col += starLayer(d, 100.0, 0.06, 17.0, 0.6);
    col += starLayer(d, 220.0, 0.025, 43.0, 0.35);
    return col * uStarBrightness;
}

// --- Matter -------------------------------------------------------------------------------

// The color of neutron star matter at log10 temperature logT and log10 density logRho.
// Films and artists show neutron stars blue-white; their true thermal glow is in X-rays
// (10^6 K and up), whose visible tail is the bluish white of any very hot blackbody. The
// dense cores are drawn white-hot, the thinner layers and tails blue; matter heated by the
// merger's shocks (10^10-10^11 K) glows brighter.
vec3 matterColor(float logT, float logRho, float doppler) {
    float heat = clamp((logT - 8.5) / 2.5, 0.0, 1.2);
    float core = smoothstep(13.6, 14.9, logRho);
    float Td = mix(12000.0, 40000.0, heat) * doppler;
    vec3 c = blackbodyHue(Td);
    if (uFilm == 1) c = mix(vec3(0.16, 0.48, 1.0), vec3(0.88, 0.95, 1.0), core);
    return c * (0.55 + 1.2 * heat + 2.2 * core);
}

// The simulation's matter at x: extinction (per render unit) and emission.
void sampleMatter(vec3 x, vec3 dir, float footprint, inout vec3 color, inout float trans, float ds) {
    vec3 uvw = (x - uBoxMin) / uBoxSize;
    if (any(lessThan(uvw, vec3(0.0))) || any(greaterThan(uvw, vec3(1.0)))) return;
    vec4 a = texture(uVolumeA, uvw);
    if (a.x < 5.0) return;
    vec4 b = texture(uVolumeB, uvw);
    vec3 label = vec3(a.z, a.w, b.w);
    // Detail finer than the grid, carried by the matter's labels: the hot surface boils
    // and the tidal tails break into strands.
    float detail = fbmFiltered(label * 7.0, 5, footprint * 7.0 * uLabelScale);
    float logRho = a.x + 0.9 * detail;
    float alpha = uOpacityScale * pow(10.0, 0.8 * (logRho - 15.0));
    if (alpha * ds < 1e-5) return;
    // Doppler factor of the flow toward the camera.
    vec3 v = b.xyz;
    float beta2 = min(dot(v, v), 0.8);
    float gamma = inversesqrt(1.0 - beta2);
    float D = 1.0 / (gamma * (1.0 - dot(v, -dir)));
    // Veins of brighter plasma along the ridges of the labels' noise.
    float veins = ridgedFiltered(label * 4.0 + 11.0, 4, footprint * 4.0 * uLabelScale);
    vec3 emit = matterColor(a.y, logRho, D) * (0.6 + 1.1 * veins * veins) * uMatterGlow * D * D * D;
    float absorbed = 1.0 - exp(-alpha * ds);
    // Thin, hot gas (tidal tails, the torus) glows more than its opacity alone would give.
    float thinGlow = 0.06 * uOpacityScale * pow(10.0, 0.45 * (logRho - 15.0)) * ds;
    color += trans * emit * (absorbed + thinGlow);
    trans *= 1.0 - absorbed;
}

// --- Gravitational waves --------------------------------------------------------------------

vec3 waveSheet(vec3 p) {
    float r = length(p.xy);
    if (r < uWaveInner || r > uWaveOuter) return vec3(0.0);
    float tc = uWaveNow - r * uWaveStep;  // retarded time, as a texture coordinate
    if (tc < 0.0) return vec3(0.0);
    vec2 H = texture(uWaves, vec2(tc, 0.5)).rg;
    float phi = atan(p.y, p.x);
    // Re[H e^(-2 i phi)]: the pattern of a rotating quadrupole.
    float w = H.x * cos(2.0 * phi) + H.y * sin(2.0 * phi);
    float crest = pow(max(w, 0.0), 5.0);
    float fadeIn = smoothstep(uWaveInner, 2.0 * uWaveInner, r);
    float fadeOut = 1.0 - smoothstep(0.5 * uWaveOuter, uWaveOuter, r);
    // 1/r amplitude, shown as 1/sqrt(r) so the outgoing rings stay visible.
    float amp = sqrt(uWaveInner / r);
    vec3 cyan = uFilm == 1 ? vec3(0.25, 0.85, 1.0) : vec3(0.55, 0.8, 1.0);
    return cyan * crest * amp * fadeIn * fadeOut * uWaveGlow;
}

// --- The debris cloud ------------------------------------------------------------------------

float wrapAngle(float a) {
    return a - 2.0 * PI * floor((a + PI) / (2.0 * PI));
}

// How the cloud's gas is spread at q = x / edge (s = |q| < 1), relative to its mean: the
// tidal tails' share and the shocked gas's.
vec2 cloudEnvelope(vec3 q, float s) {
    // The tidal tails: flung from the stars' far sides as they touched and wound back by
    // the turning binary into two trailing arms (the orbit turns counterclockwise about z),
    // thin near the orbital plane and fanning out as they go.
    float phi = atan(q.y, q.x);
    float arm = uCloudMore.y - 5.0 * s;
    float d0 = wrapAngle(phi - arm);
    float d1 = wrapAngle(phi - arm - PI);
    float w = 0.3 + 0.55 * s;
    float arms = exp(-d0 * d0 / (2.0 * w * w)) + exp(-d1 * d1 / (2.0 * w * w));
    float h = 0.06 + 0.2 * s;
    float tails = arms * exp(-q.z * q.z / (2.0 * h * h));
    // The shocked gas squeezed out from between the cores: all round, more toward the poles.
    float c = abs(q.z) / max(s, 1e-4);
    float shock = 0.18 + 0.7 * c * c;
    return vec2(4.0 * tails, shock);
}

vec3 cloudColor(float T) {
    vec3 c = blackbodyHue(T);
    if (uFilm == 1) {
        // Films paint hot debris as fire: white at its heart, then gold, orange, and a dull
        // red where it is coolest.
        float h = clamp((T - 1500.0) / 9000.0, 0.0, 1.0);
        vec3 fire = mix(vec3(0.6, 0.09, 0.03), vec3(1.0, 0.42, 0.08), smoothstep(0.0, 0.35, h));
        fire = mix(fire, vec3(1.0, 0.76, 0.36), smoothstep(0.3, 0.7, h));
        fire = mix(fire, vec3(1.0, 0.95, 0.86), smoothstep(0.7, 1.0, h));
        c = mix(c, fire, 0.75);
    }
    return c;
}

// The baked noise at p (in noise cells, four to the texture), filtered over a footprint of
// fp cells.
vec4 cloudNoise(vec3 p, float fp) {
    return textureLod(uNoise, p * 0.25, log2(max(fp * 32.0, 1.0)));
}

// The billows' density (about 0-2) at p (noise cells), from the baked noise there.
float lumpsAt(vec4 n) {
    return 0.45 + 0.9 * n.r + 0.6 * n.g;
}

void sampleCloud(vec3 x, float footprint, inout vec3 color, inout float trans, float ds) {
    float R = uCloud.x;
    float r = length(x);
    if (r >= R || r < 0.5 * uCloudMore.x) return;
    float s = r / R;
    vec3 q = x / R;
    vec3 n = x / max(r, 1e-6);
    // The ragged hole the remnant sits in.
    float inner = smoothstep(0.5 * uCloudMore.x, 1.6 * uCloudMore.x, r);
    vec2 parts = cloudEnvelope(q, s) * inner;
    float envelope = parts.x + parts.y;
    if (envelope < 0.02) return;
    // Billows: turbulence (Kelvin-Helmholtz between the tails and the shocked gas, Rayleigh-
    // Taylor as the hot gas pushes out) frozen into the coasting gas, so it grows as the cloud
    // does, and still churning as it expands; finer curls on the big ones.
    float fp = footprint / R * 3.4;
    vec3 p = q * 3.4 + vec3(0.0, 0.0, uCloudMore.z);
    vec4 w = cloudNoise(p * 0.5 + 1.7, fp * 0.5);
    p += 1.1 * vec3(w.b, w.a, w.r);
    vec4 big = cloudNoise(p, fp);
    float fine = cloudNoise(p * 4.3 + 7.3, fp * 4.3).r;
    float lumps = lumpsAt(big) + 0.35 * fine;
    // The gas gathers into billows with clear gaps between them, more of them where the
    // envelope is dense (the tails), fewer as the cloud thins. Their mean optical depth stays
    // the cloud's: the fewer they are, the denser each.
    float fill = clamp(0.3 + 0.12 * log(max(uCloud.y, 1e-3)) + 0.25 * min(envelope, 1.5), 0.08, 0.8);
    float threshold = 1.15 - 0.95 * fill;
    float edge = 1.0 - smoothstep(0.5 + 0.35 * big.r, 1.0, s);
    float density = envelope * edge * smoothstep(threshold, threshold + 0.18, lumps) / fill;
    if (density < 1e-3) return;
    float alpha = uCloud.y * density / R;
    // Lit from within by the hot remnant: its light dims through the cloud on the way out,
    // and each billow is shaded by the one just inside it, so the cloud glows through its
    // gaps and its billows stand out dark-edged, as fire seen through its own smoke.
    float shade = exp(-2.5 * max(lumpsAt(cloudNoise(p - n * 0.3, fp)) - threshold, 0.0));
    float reach = exp(-0.45 * min(uCloud.y, 6.0) * s);
    float light = 0.22 + (0.4 + 1.6 * uCloudMore.w * exp(-3.0 * s)) * shade * reach;
    // Hotter deeper in and on the billows' crests. The tidal tails were flung off, not
    // shocked: they stay cooler (and so redder) than the shocked gas.
    float shocked = parts.y / envelope;
    float T = uCloud.w * (0.75 + 0.5 * (1.0 - s)) * (0.85 + 0.3 * lumps) * mix(0.55, 1.15, shocked);
    // Thinning gas still glows with the heat of the decays inside it: the glow per optical
    // depth rises as the optical depth falls, so the last wisps stay visible.
    float thin = pow(max(1.0 / max(uCloud.y, 1e-3), 1.0), 0.5);
    vec3 emit = cloudColor(T) * uCloud.z * (0.3 + 0.9 * lumps * lumps) * light * thin;
    float absorbed = 1.0 - exp(-alpha * ds);
    color += trans * emit * absorbed;
    trans *= 1.0 - absorbed;
}

// --- Kilonova --------------------------------------------------------------------------------

// The share of a component's mass in direction n (unit): blue toward the poles (driven off
// the hot remnant by shocks and neutrinos), purple everywhere, red in the orbital plane
// (the tidal tails), with the simulation's ejecta pattern on top.
// How far out component k reaches in direction n, in units of its edge radius: the blue
// stretched toward the poles, the red flattened into the orbital plane (src/merger.h).
float ejectaStretch(int k, vec3 n) {
    vec2 ab = k == 0 ? vec2(0.8, 1.15) : (k == 1 ? vec2(1.15, 0.85) : vec2(1.6, 0.55));
    float c2 = n.z * n.z;
    return inversesqrt((1.0 - c2) / (ab.x * ab.x) + c2 / (ab.y * ab.y));
}

float ejectaAngular(int k, vec3 n) {
    float c = abs(n.z);
    if (k == 0) return 0.06 + 1.8 * smoothstep(0.3, 0.85, c);
    if (k == 1) return 0.45 + 0.55 * (1.0 - c * c);
    float th = acos(clamp(n.z, -1.0, 1.0));
    float ph = atan(n.y, n.x) + PI;
    float sim = texture(uEjectaSky, vec2(ph / (2.0 * PI), th / PI)).r;
    return (0.2 + 2.2 * exp(-c * c / (0.3 * 0.3))) * (0.5 + 2.5 * sim);
}

vec3 kilonovaColor(int k, float T) {
    vec3 c = blackbodyHue(T);
    if (uFilm == 1) {
        // Films paint the lanthanide-poor ejecta electric blue and the lanthanide-rich ones
        // fiery orange.
        vec3 film = k == 0 ? vec3(0.22, 0.58, 1.0) : (k == 1 ? vec3(1.0, 0.42, 0.42) : vec3(1.0, 0.45, 0.1));
        c = mix(c, film, 0.65);
    }
    return c;
}

void sampleKilonova(vec3 x, float footprint, inout vec3 color, inout float trans, float ds) {
    // All three components coast out together from the same place, so their structure is one
    // pattern frozen in velocity space (x / t): measured here against the fastest one's edge.
    // Rayleigh-Taylor plumes (the r-process heat pushing through the ejecta) give it billows
    // like a cauliflower's.
    float Rmax = max(uKilonova[0].x * 1.15, max(uKilonova[1].x * 1.15, uKilonova[2].x * 1.6));
    float r = length(x);
    if (Rmax <= 0.0 || r >= Rmax) return;
    vec3 q = x / Rmax;
    float fp = footprint / Rmax;
    float billows = cloudNoise(q * 4.0, fp * 4.0).r;
    float plumes = cloudNoise(q * 9.0 + 7.3, fp * 9.0).g;
    float lumps = 0.5 + 0.9 * billows + 0.5 * plumes;
    vec3 n = x / max(r, 1e-6);
    for (int k = 0; k < 3; ++k) {
        vec4 kn = uKilonova[k];
        if (kn.x <= 0.0) continue;
        float s = r / (kn.x * ejectaStretch(k, n));
        if (s >= 1.0) continue;
        // The edge is where the billows give out.
        float edge = smoothstep(1.0, 0.55 + 0.3 * billows, s);
        float shape = edge * (1.0 - 0.6 * s * s) * ejectaAngular(k, n);
        float density = shape * smoothstep(0.15, 0.85, lumps);
        // The real ejecta are opaque for days: shown at no more than ~25 optical depths so the
        // billowing surface, not just a sharp edge, can be seen.
        float alpha = min(kn.w, 25.0) * 2.2 * density / kn.x;
        if (alpha * ds < 1e-6) continue;
        // Brighter on the plumes' crests and deeper in (the heat comes from the decays inside).
        float glow = (0.35 + 1.4 * lumps * lumps) * (0.7 + 0.6 * (1.0 - s));
        // Once the ejecta thin out they still give off all the decays' heat: the glow per unit
        // optical depth rises as the optical depth falls, so a thin nebula stays visible.
        float thin = pow(max(1.0 / max(kn.w * 2.2, 1e-3), 1.0), 0.75);
        vec3 emit = kilonovaColor(k, kn.y * (0.85 + 0.3 * lumps)) * kn.z * glow * thin;
        float absorbed = 1.0 - exp(-alpha * ds);
        color += trans * emit * absorbed;
        trans *= 1.0 - absorbed;
    }
}

// --- Jet and cocoon ------------------------------------------------------------------------------

void sampleJet(vec3 x, float footprint, inout vec3 color, inout float trans, float ds) {
    if (uJet.z <= 0.0) return;
    float along = dot(x, uJetAxis);
    float side = sign(along);
    float z = abs(along);
    vec3 perp = x - along * uJetAxis;
    float rho = length(perp);
    // The jet: a narrow cone from the hole to its head, with knots moving out.
    if (z < uJet.x && z > 0.0) {
        float width = max(z * tan(uJet.w), 0.002 * uJet.x);
        float core = exp(-(rho * rho) / (width * width));
        if (core > 1e-3) {
            float knots = 0.6 + 0.8 * fbmFiltered(vec3(rho / width, z / uJet.x * 12.0 - uTime * 3.0, side * 5.0), 3, footprint / width);
            float head = 1.0 + 3.0 * exp(-pow((uJet.x - z) / (0.05 * uJet.x), 2.0));
            vec3 c = uFilm == 1 ? vec3(1.0, 0.55, 0.18) : vec3(0.75, 0.85, 1.0);
            float a = core * knots * head * uJet.z / max(width, 1e-6) * 0.5;
            color += trans * c * a * ds;
        }
    }
    // The cocoon: after breakout, hot gas spills sideways from where the jet left the
    // ejecta, a bubble swelling around each breakout point, brightest at its shocked rim.
    if (uJet.y > 0.0) {
        vec3 center = side * uJetMore.x * uJetAxis;
        float d = length(x - center);
        float shell = exp(-pow((d - uJet.y) / (0.12 * uJet.y), 2.0));
        float inner = smoothstep(uJet.y, 0.3 * uJet.y, d) * 0.15;
        if (shell + inner < 1e-4) return;
        float rough = 0.5 + cloudNoise((x - center) / uJet.y * 4.0, footprint / uJet.y * 4.0).r;
        float alpha = (shell + inner) * rough * 0.9 / uJet.y;
        // The cocoon shines in X-rays as it cools; in film colors it is drawn violet, as
        // X-ray images are.
        vec3 c = uFilm == 1 ? vec3(0.62, 0.35, 1.0) : vec3(0.7, 0.8, 1.0);
        float absorbed = 1.0 - exp(-alpha * ds * 0.6);
        color += trans * c * uJetMore.y * absorbed;
        trans *= 1.0 - absorbed * 0.5;
    }
}

// --- Main ----------------------------------------------------------------------------------------

void main() {
    vec2 ndc = ((gl_FragCoord.xy - uViewport) / uResolution) * 2.0 - 1.0;
    ndc.x *= uResolution.x / uResolution.y;
    vec3 dir = normalize(uCamBasis * vec3(ndc * uTanHalfFov, 1.0));
    vec3 x = uCamPos;

    // How far out anything is: past this, rays only bend a little more on their way out.
    float reach = 0.0;
    if (uVolumeOn == 1) reach = max(reach, length(uBoxMin + 0.5 * uBoxSize) + 0.9 * uBoxSize);
    for (int k = 0; k < 3; ++k) reach = max(reach, uKilonova[k].x * 1.6);
    reach = max(reach, uCloud.x);
    if (uJet.z > 0.0) reach = max(reach, max(uJet.x, uJetMore.x + 1.3 * uJet.y));
    if (uWavesOn == 1) reach = max(reach, uWaveOuter);
    for (int k = 0; k < uLensCount; ++k) reach = max(reach, length(uLens[k].xyz) + 30.0 * uLens[k].w);

    vec3 color = vec3(0.0);
    float trans = 1.0;
    bool captured = false;
    float cell = uVolumeOn == 1 ? uBoxSize / 128.0 : 1e30;
    float travelled = 0.0;
    for (int i = 0; i < 900; ++i) {
        float r = length(x);
        if (r > reach * 1.02 && dot(x, dir) > 0.0) break;
        // Jump straight to the scene from far away.
        if (r > reach * 1.05) {
            float b = dot(x, dir);
            float c = dot(x, x) - reach * reach;
            float disc = b * b - c;
            if (disc < 0.0) break;
            float t = -b - sqrt(disc);
            x += dir * max(t, 0.0);
            travelled += max(t, 0.0);
            continue;
        }
        // Step: small near the lenses and inside the matter.
        float ds = max(0.03 * reach, 1e-4);
        for (int k = 0; k < uLensCount; ++k) {
            float rk = length(x - uLens[k].xyz);
            ds = min(ds, max(0.12 * rk, 0.02 * uLens[k].w));
        }
        vec3 uvw = (x - uBoxMin) / uBoxSize;
        bool inBox = uVolumeOn == 1 && all(greaterThan(uvw, vec3(-0.02))) && all(lessThan(uvw, vec3(1.02)));
        if (inBox) ds = min(ds, 0.8 * cell);
        for (int k = 0; k < 3; ++k)
            if (uKilonova[k].x > 0.0 && r < uKilonova[k].x * 1.62) ds = min(ds, 0.012 * uKilonova[k].x);
        if (uCloud.x > 0.0 && r < uCloud.x * 1.02) ds = min(ds, 0.013 * uCloud.x);
        if (uJet.z > 0.0 && r < reach) ds = min(ds, 0.01 * max(uJet.x, 1e-4) + 0.05 * abs(length(x - dot(x, uJetAxis) * uJetAxis)));

        // Bend: midpoint step of x'' = -3 M h^2 x / r^5 for each lens.
        vec3 mid = x + 0.5 * ds * dir;
        vec3 bend = vec3(0.0);
        for (int k = 0; k < uLensCount; ++k) {
            vec3 rel = mid - uLens[k].xyz;
            float rk = length(rel);
            vec3 hv = cross(rel, dir);
            bend -= 3.0 * uLens[k].w * dot(hv, hv) * rel / pow(rk, 5.0);
            if (k == uHole && rk < 2.0 * uLens[k].w) captured = true;
        }
        if (captured) break;
        vec3 newDir = normalize(dir + bend * ds);
        vec3 newX = x + 0.5 * ds * (dir + newDir);

        float footprint = uPixelAngle * (travelled + 0.5 * ds);
        if (uWavesOn == 1 && sign(newX.z) != sign(x.z)) {
            float f = x.z / (x.z - newX.z);
            vec3 p = mix(x, newX, f);
            color += trans * waveSheet(p) * min(1.0 / max(abs(dir.z), 0.08), 6.0) * 0.25;
        }
        vec3 m = 0.5 * (x + newX);
        if (inBox) sampleMatter(m, dir, footprint, color, trans, ds);
        if (uCloud.x > 0.0) sampleCloud(m, footprint, color, trans, ds);
        if (uKilonovaOn == 1) sampleKilonova(m, footprint, color, trans, ds);
        sampleJet(m, footprint, color, trans, ds);

        x = newX;
        dir = newDir;
        travelled += ds;
        if (trans < 0.004) {
            trans = 0.0;
            break;
        }
    }
    if (!captured && trans > 0.0) color += trans * sky(dir);
    vec3 toCenter = normalize(-uCamPos);
    float ang = acos(clamp(dot(normalize(uCamBasis * vec3(ndc * uTanHalfFov, 1.0)), toCenter), -1.0, 1.0));
    // The gamma-ray flash, seen as a glare around the jet's base.
    if (uFlash > 0.0) color += vec3(1.0, 0.85, 0.7) * uFlash * (0.02 / (ang * ang + 0.0004) * 0.002 + 0.25 * exp(-ang * 6.0));
    // The collision's flash: a white-hot glare the size of the merging cores, and a wide halo.
    if (uMergeFlash > 0.0) {
        float core = max(1.5 * uCloudMore.x / max(length(uCamPos), 1e-6), 0.03);
        color += vec3(1.0, 0.95, 0.88) * uMergeFlash * (1.2 * exp(-ang * ang / (2.0 * core * core)) + 0.3 * exp(-ang * 4.0));
    }
    outColor = vec4(color, 1.0);
}
