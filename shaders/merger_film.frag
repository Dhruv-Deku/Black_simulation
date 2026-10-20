#version 460 core

// The neutron star merger as the reference film shows it (src/film.h): its shots one after
// another, each matched to the reference frame by frame: where things are, how big, what
// colour, when. Positions and sizes are written in the reference's own pixels (380 x 688,
// y down) through ref() and px(), so they can be checked against it directly.
//
// Output is linear light; the post pass adds bloom, clamps and encodes to sRGB, so a colour
// read off the reference as sRGB is written srgb(r, g, b), and white-hot parts go above 1.

in vec2 vUV;
layout(location = 0) out vec4 outColor;

uniform vec2 uResolution;    // the target's size in pixels
uniform float uFrameHeight;  // the 9:16 frame's height in pixels (it is centred)
uniform float uTime;         // seconds into the film, 0 to 17.6
// The binary of the first and last shots (src/film.cpp), camera space.
uniform vec3 uStar0, uStar1;
uniform float uGlare;        // how far each star's white glare reaches
uniform float uMerge;        // 0 apart, 1 fused
uniform vec3 uOrbitCentre;
uniform float uOrbitPhase, uOrbitRadius;

#include "noise.glsl"

const float PI = 3.14159265358979;
const float FOCAL = 1.2071068;  // 45-degree vertical field of view, frame heights
const vec3 ORBIT_U = vec3(0.0, 1.0, 0.0);
const vec3 ORBIT_W = vec3(-0.42, 0.0, 0.9075241);
const vec3 ORBIT_N = vec3(0.9075241, 0.0, 0.42);  // cross(U, W)
const float ORBIT_COS = 0.42;

// When each shot starts (s).
const float T_FLASH = 6.63, T_WARP = 6.69, T_STREAK = 6.86, T_BLUR = 6.99, T_WAVES = 7.15;
const float T_BURST = 8.07, T_MERGE = 9.02, T_WHIP = 9.97, T_HEAVY = 10.07, T_FADE = 10.9;
const float T_SHELL = 11.25, T_BLACK = 14.72, T_CALM = 14.9;

// The axis the jets take, upper left to lower right.
const vec2 JET = vec2(-0.848, 0.530);

float pixel;  // one output pixel, in frame heights

// The reference's pixel (x, y) and a length of n of its pixels, in frame coordinates.
vec2 ref(float x, float y) { return vec2(x - 190.0, 344.0 - y) / 688.0; }
float px(float n) { return n / 688.0; }

vec3 srgb(float r, float g, float b) { return pow(vec3(r, g, b) / 255.0, vec3(2.2)); }

float saturate(float x) { return clamp(x, 0.0, 1.0); }
float band(float a, float b, float t) { return smoothstep(a, b, t); }
float sq(float x) { return x * x; }
float gauss(float x) { return exp(-x * x); }

// --- Shared pieces ---------------------------------------------------------------------------

// The faint stars behind everything: one in about every eighth cell of a fine grid, a few
// bright, most barely there.
vec3 starfield(vec2 p, float amount) {
    vec3 col = vec3(0.0);
    for (int layer = 0; layer < 2; ++layer) {
        float cell = layer == 0 ? px(9.0) : px(16.0);
        vec2 q = p / cell;
        vec2 base = floor(q);
        for (int y = -1; y <= 1; ++y)
            for (int x = -1; x <= 1; ++x) {
                vec2 c = base + vec2(x, y);
                vec4 h = hash43(vec3(c, 3.7 + float(layer) * 11.0));
                if (h.w > (layer == 0 ? 0.13 : 0.09)) continue;
                vec2 at = (c + 0.2 + 0.6 * h.xy) * cell;
                float d = length(p - at);
                float size = max(pixel * 0.55, px(0.45));
                float b = (0.15 + 1.6 * pow(h.z, 6.0)) * (layer == 0 ? 0.7 : 1.0);
                col += vec3(0.85, 0.9, 1.0) * b * exp(-d * d / (2.0 * size * size));
            }
    }
    return col * amount;
}

// Distance from p to the segment a-b, and how far along it (0 to 1) p lies.
float segment(vec2 p, vec2 a, vec2 b, out float along) {
    vec2 ab = b - a;
    along = clamp(dot(p - a, ab) / dot(ab, ab), 0.0, 1.0);
    return length(p - a - ab * along);
}

// Distance to an ellipse's outline (semi-axes ra along u, rb across), approximately but with
// an even line width all round: the implicit function over its gradient.
float ellipseLine(vec2 p, vec2 centre, vec2 u, float ra, float rb) {
    vec2 d = p - centre;
    vec2 q = vec2(dot(d, u), dot(d, vec2(-u.y, u.x)));
    float k = length(q / vec2(ra, rb));
    float g = length(q / vec2(ra * ra, rb * rb)) / max(k, 1e-6);
    return (k - 1.0) / max(g, 1e-6);
}

vec3 viewRay(vec2 p) { return normalize(vec3(p, -FOCAL)); }

// --- Billowing clouds: the kilonova's ejecta, the red remnant -----------------------------------

// Billows: the folded noise |n| is rounded on top and creased where n crosses zero, as
// cumulus are; summed over three scales.
float billow(vec3 q) {
    float s = abs(gradientNoise(q)) * 0.5;
    s += abs(gradientNoise(q * 2.2 + 7.1)) * 0.32;
    s += abs(gradientNoise(q * 4.9 + 3.7)) * 0.18;
    return s;  // about 0 (creases) to 0.8 (puff tops)
}

struct Cloud {
    vec3 centre;     // camera space
    float radius;
    vec3 axis;       // flattened along this axis by `squash`
    float squash;
    float spin;      // its puffs turn about the vertical by this angle
    float seed;
    float rough;     // how far its puffs stand out, in radii
    float scale;     // puffs across a radian
};

vec3 cloudLocal(vec3 x, Cloud c) {
    vec3 d = x - c.centre;
    d += c.axis * dot(d, c.axis) * (1.0 / c.squash - 1.0);
    float s = sin(c.spin), co = cos(c.spin);
    return vec3(co * d.x + s * d.z, d.y, -s * d.x + co * d.z);
}

// The puffs are a function of direction only, so they stand out radially, as an explosion's
// debris does; lift is in radii.
float cloudLift(vec3 n, Cloud c) {
    return c.rough * (billow(n * c.scale + c.seed) - 0.33);
}

float cloudDistance(vec3 x, Cloud c, out float lift) {
    vec3 d = cloudLocal(x, c);
    float len = length(d);
    lift = cloudLift(d / max(len, 1e-6), c);
    return (len - c.radius * (1.0 + lift)) * c.squash * 0.5;
}

struct CloudLook {
    vec3 deep, mid, lit, rim;
    float grain;     // fine speckle in the creases
};

// The cloud along the ray rd from the camera: colour and coverage. It glows by itself (the
// ejecta are hot): bright where it puffs out, deep blue in the creases, a bright rim.
vec4 drawCloud(vec3 rd, Cloud c, CloudLook look) {
    float reach = c.radius * (1.0 + c.rough * 0.5) * 1.05;
    float b = dot(rd, c.centre);
    float disc = b * b - dot(c.centre, c.centre) + reach * reach;
    if (disc <= 0.0) return vec4(0.0);
    float sq_ = sqrt(disc);
    float t = max(b - sq_, 0.0), tEnd = b + sq_;
    float closest = 1e9;
    float lift = 0.0;
    bool hit = false;
    for (int i = 0; i < 72; ++i) {
        float d = cloudDistance(rd * t, c, lift);
        closest = min(closest, d / t);
        if (d < 0.0002 * t) {
            hit = true;
            break;
        }
        t += max(d, 0.0015 * c.radius);
        if (t > tEnd) break;
    }
    // A soft, ragged edge: coverage from how near the ray came, in pixels.
    float cover = hit ? 1.0 : saturate(1.0 - closest / (pixel / FOCAL * 1.5));
    if (cover <= 0.0) return vec4(0.0);
    vec3 x = rd * t;
    float l2;
    float e = 0.005 * c.radius;
    vec2 k = vec2(1.0, -1.0);
    vec3 n = normalize(k.xyy * cloudDistance(x + k.xyy * e, c, l2) + k.yyx * cloudDistance(x + k.yyx * e, c, l2) +
                       k.yxy * cloudDistance(x + k.yxy * e, c, l2) + k.xxx * cloudDistance(x + k.xxx * e, c, l2));
    cloudDistance(x, c, lift);
    // Creases (low lift) deep, puffs bright; lit a little from the upper left.
    float height = saturate((lift / max(c.rough, 1e-3) + 0.33) / 0.7);
    float ao = pow(height, 0.7);
    vec3 L = normalize(vec3(-0.5, 0.6, 0.62));
    float diff = max(dot(n, L), 0.0);
    float facing = max(dot(n, -rd), 0.0);
    float rim = sq(1.0 - facing);
    float speck = valueNoise(cloudLocal(x, c) / c.radius * 70.0);
    float grain = 1.0 + look.grain * (speck - 0.5) * (1.5 - height);
    vec3 col = mix(look.deep, look.mid, ao) * (0.8 + 0.4 * diff) + look.lit * pow(ao, 2.2) * (0.55 + 0.65 * diff);
    col = col * grain + look.rim * rim;
    return vec4(col, cover);
}

// A cloud drawn at frame position f with radius r (frame heights), at the camera's distance 1.
Cloud cloudAt(vec2 f, float r) {
    Cloud c;
    c.centre = vec3(f / FOCAL, -1.0);
    c.radius = r / FOCAL;
    c.axis = vec3(0.0, 1.0, 0.0);
    c.squash = 1.0;
    c.spin = 0.0;
    c.seed = 0.0;
    c.rough = 0.3;
    c.scale = 5.0;
    return c;
}

// --- Shot 1: two neutron stars spiral together (0 - 6.63 s) ------------------------------------

// The wakes the stars trail, as in the reference: a crescent sheet of glowing gas hugging each
// star's outer side, reaching far back along the orbit behind it (trailing, spiralling out)
// and a little ahead; broad and bright while the star is on the near side of the orbit, a
// thin bright arc while it is on the far side. Older, thinner arcs lie further out.
vec3 wakes(vec2 p, float t) {
    vec3 rd = viewRay(p);
    float denom = dot(rd, ORBIT_N);
    if (abs(denom) < 1e-4) return vec3(0.0);
    float s = dot(uOrbitCentre, ORBIT_N) / denom;
    if (s <= 0.0) return vec3(0.0);
    vec3 q = rd * s - uOrbitCentre;
    vec2 uw = vec2(dot(q, ORBIT_U), dot(q, ORBIT_W)) / uOrbitRadius;
    float rho = length(uw);
    float phi = atan(uw.y, uw.x);
    // How long a step in rho looks on the screen here (the orbit is drawn squeezed sideways).
    float onScreen = sqrt(sq(ORBIT_COS * sin(phi)) + sq(cos(phi)));
    float starR = uGlare / uOrbitRadius;  // the star's radius in orbit radii
    float grow = band(0.0, 5.8, t);
    vec3 col = vec3(0.0);
    for (int arm = 0; arm < 2; ++arm) {
        float psi = uOrbitPhase + PI * float(arm);
        float nearness = 0.5 + 0.5 * sin(psi);  // of the star itself
        float delta = mod(phi - psi + PI, 2.0 * PI) - PI;  // < 0: behind the star
        float behind = max(-delta, 0.0);
        // The sheet's screen width and where it lies: just outside the star's limb.
        float width = mix(0.06, 0.3, nearness) * (1.0 + 0.8 * behind) * (1.0 + 0.6 * grow);
        float rc = 1.0 + (0.9 * starR + 0.5 * width) / onScreen + 0.14 * behind;
        float x = (rho - rc) * onScreen / (0.5 * width);
        float reach = delta < 0.0 ? exp(-sq(delta / (1.2 + 0.6 * grow))) : exp(-sq(delta / 0.4));
        float body = gauss(x) * (0.6 + 0.4 * band(0.8, -0.8, x));
        float edge = gauss((x + 0.8) / 0.22);
        float streak = 0.55 + 0.9 * fbm(vec3(x * 1.3, delta * 2.0 - t * 0.15, float(arm) * 5.0));
        float bright = mix(0.6, 1.0, nearness);
        col += (srgb(78, 188, 214) * body * streak + srgb(175, 238, 245) * edge * 0.55) * reach * bright;
        // An older, thinner arc further out, longer.
        float rc2 = rc + 0.45 + 0.2 * behind;
        float x2 = (rho - rc2) * onScreen / px(7.0) * uOrbitRadius * FOCAL;
        float reach2 = delta < 0.0 ? exp(-sq(delta / 1.4)) : exp(-sq(delta / 0.5));
        col += srgb(70, 160, 195) * gauss(x2) * reach2 * (0.04 + 0.4 * grow);
    }
    return col * (0.6 + 1.6 * pow(grow, 1.2));
}

// The stars: a white disk that the camera's sensor saturates, a tight white-cyan glare and a
// wide teal glow. Drawn as one shape, so that as they touch they fuse into a peanut.
vec3 binaryStars(vec2 p, vec3 bg, float glow) {
    vec2 c0 = uStar0.xy * (FOCAL / -uStar0.z), c1 = uStar1.xy * (FOCAL / -uStar1.z);
    float r0 = uGlare * FOCAL / -uStar0.z, r1 = uGlare * FOCAL / -uStar1.z;
    float d0 = length(p - c0) - 0.97 * r0;
    float d1 = length(p - c1) - 0.97 * r1;
    // Smooth union: the waist fills in a little as they merge.
    float k = mix(0.02, 0.15, uMerge) * max(r0, r1);
    float hh = saturate(0.5 + 0.5 * (d1 - d0) / k);
    float sd = mix(d1, d0, hh) - k * hh * (1.0 - hh);
    float r = max(r0, r1);
    float inner = exp(-max(sd, 0.0) / (0.22 * r));
    float outer = exp(-max(sd, 0.0) / (1.3 * r));
    vec3 col = bg + glow * (srgb(165, 225, 245) * 1.0 * inner + srgb(35, 125, 175) * 0.4 * outer);
    float diskMask = smoothstep(pixel, -pixel, sd);
    vec3 disk = mix(srgb(215, 240, 255) * 1.3, vec3(3.0), saturate(-sd / (0.2 * r)));
    return mix(col, disk, diskMask);
}

vec3 shotInspiral(vec2 p, float t) {
    // The dark navy sky fills with teal light, greying as they merge.
    float fog = 0.75 * pow(band(2.0, 6.2, t), 1.5);
    vec2 c = ref(176.0, 320.0);
    float vignette = exp(-dot(p - c, p - c) / 0.25);
    vec3 fogCol = mix(srgb(16, 62, 105), srgb(30, 92, 118), band(5.6, 6.3, t));
    vec3 bg = srgb(2, 5, 12) + fog * (fogCol * (0.6 + 0.4 * vignette) + srgb(10, 30, 45) * fbm(vec3(p * 3.0, t * 0.2)));
    bg += starfield(p, 1.0 - 0.7 * fog);
    bg += wakes(p, t);
    // Grey smoke torn off as they merge, with darker gaps.
    float smoke = band(5.6, 6.1, t);
    float n = fbm(vec3(p.x * 4.0, p.y * 1.5, t * 0.6));
    bg = mix(bg, bg * 0.55, smoke * band(0.55, 0.75, n) * 0.8);
    bg += srgb(70, 105, 120) * smoke * 0.5 * n * exp(-length(p - c) / 0.35);
    // The flash: the glare swells over its last tenth of a second.
    float flash = band(6.5, 6.63, t);
    vec3 col = binaryStars(p, bg, 1.0 + 5.0 * flash * flash);
    col += flash * flash * srgb(20, 160, 175) * 1.5;
    return col;
}

// --- The flash and the cuts that follow (6.63 - 7.15 s) -----------------------------------------

vec3 shotFlash(vec2 p, float t) {
    vec2 c = ref(180.0, 320.0);
    float d = length(p - c);
    // Brightest along the jets' diagonal, as in the reference.
    vec2 q = p - c;
    float diag = gauss(dot(q, vec2(-JET.y, JET.x)) / 0.22);
    vec3 col = srgb(4, 62, 85) + srgb(6, 110, 115) * diag + srgb(10, 190, 220) * 1.3 * exp(-d / 0.17);
    col += vec3(3.0) * smoothstep(px(56.0), px(50.0), d) + srgb(120, 255, 255) * 1.3 * exp(-max(d - px(52.0), 0.0) / px(18.0));
    col += srgb(8, 50, 60) * fbm(vec3(p.x * 6.0, p.y * 1.2, 3.0));
    return col;
}

// Gravitational-wave rings: circles in the orbital plane, perpendicular to the jets, seen
// nearly edge-on (a third as wide as long), each a soft band.
vec3 waveRings(vec2 p, vec2 centre, vec2 jet, float first, float spacing, int count, vec3 tint, float width) {
    vec3 col = vec3(0.0);
    vec2 across = vec2(-jet.y, jet.x);
    for (int n = 0; n < count; ++n) {
        float rr = first + spacing * float(n);
        if (rr < px(20.0)) continue;
        float d = ellipseLine(p, centre, across, rr, rr * 0.32);
        vec2 q = p - centre;
        float side = 0.75 + 0.25 * sin(atan(dot(q, jet), dot(q, across)) + 0.6);
        float fade = band(px(15.0), px(60.0), rr) * exp(-rr / px(500.0));
        col += tint * side * fade * (gauss(d / width) * 0.8 + gauss(d / (width * 3.0)) * 0.35);
    }
    return col;
}

// Short white streaks flying out from c, as debris rushing past the camera.
vec3 speedLines(vec2 p, vec2 c, float t, float amount) {
    vec2 q = p - c;
    float d = length(q);
    float ang = atan(q.y, q.x);
    float bins = 110.0;
    float bin = floor(ang / (2.0 * PI) * bins);
    vec3 col = vec3(0.0);
    for (int k = 0; k < 3; ++k) {
        vec4 h = hash43(vec3(bin, 5.0 + float(k), 1.0));
        if (h.y > 0.5) continue;
        float centreAng = (bin + 0.2 + 0.6 * h.z) / bins * 2.0 * PI;
        float across = d * abs(sin(ang - centreAng));
        float rad = mix(0.06, 0.7, h.x) + (t - T_WARP) * 1.2;
        float len = px(3.0 + 5.0 * h.w);
        col += vec3(0.7, 0.85, 0.92) * gauss(across / (pixel * 0.7)) * gauss((d - rad) / len) * (0.3 + 0.5 * h.w);
    }
    return col * amount;
}

// Fluffy puffs of gas blown out both ways along the jet axis from c, out to `reach`.
float puffs(vec2 p, vec2 c, float reach, float width, float seed) {
    vec2 q = p - c;
    float a = dot(q, JET), x = dot(q, vec2(-JET.y, JET.x));
    float len = abs(a);
    float lumpy = fbm(vec3(q / width * 1.3, seed + sign(a)));
    float w = width * (0.6 + 0.8 * lumpy) * (1.0 - 0.4 * band(0.0, reach, len));
    return smoothstep(w, w * 0.3, abs(x)) * band(width * 0.4, width * 1.2, len) * (1.0 - band(reach * 0.75, reach, len)) *
           (0.5 + 0.9 * fbm(vec3(q / width * 3.0, seed + 5.0)));
}

vec3 shotWarp(vec2 p, float t) {
    vec2 c = ref(181.0, 320.0);
    vec2 q = p - c;
    float d = length(q);
    float diag = gauss(dot(q, vec2(-JET.y, JET.x)) / 0.3);
    vec3 col = srgb(6, 19, 38) + srgb(10, 30, 40) * diag + srgb(8, 22, 30) * exp(-d / 0.15);
    // Broad blurred arcs of the waves sweeping past, mostly on the right.
    col += waveRings(p, c, JET, px(250.0), px(120.0), 2, srgb(28, 80, 120), px(18.0)) * (0.45 + 0.55 * band(-0.1, 0.2, q.x));
    col += speedLines(p, c, t, 0.6);
    // Small cyan jets, puffy.
    col += srgb(70, 205, 228) * puffs(p, c, px(70.0), px(21.0), 2.0) * 1.3;
    // The new-born remnant's spark.
    col += vec3(3.0) * smoothstep(px(13.0), px(10.0), d) + srgb(60, 235, 245) * 1.8 * exp(-d / px(20.0)) +
           srgb(20, 140, 160) * 0.6 * exp(-d / px(70.0));
    return col;
}

vec3 shotStreak(vec2 p, float t) {
    vec2 c = ref(183.0, 316.0);
    vec3 col = srgb(5, 7, 11) + srgb(6, 10, 16) * exp(-dot(p - c, p - c) / 0.05);
    col += waveRings(p, c, JET, px(215.0), px(110.0), 3, srgb(22, 45, 70), px(16.0));
    // A fibrous, clumpy jet, yellow-green.
    float along;
    vec2 a0 = ref(92.0, 262.0), a1 = ref(275.0, 368.0);
    float wob = (fbm(vec3(p * 22.0, 4.0)) - 0.5) * px(16.0);
    float d = segment(p, a0, a1, along) + wob;
    float width = px(6.0) + px(8.0) * abs(along - 0.5) * 2.0;
    float fib = fbm(vec3(along * 25.0, d / px(3.0), 4.0)) * fbm(vec3(p * 60.0, 1.0)) * 1.6;
    float clumps = 0.4 + 1.2 * fbm(vec3(along * 9.0, 2.0, 7.0));
    col += srgb(160, 165, 60) * gauss(d / width) * (0.25 + 1.4 * fib) * clumps * band(0.0, 0.08, along) * band(1.0, 0.8, along);
    float dc = length(p - c);
    col += vec3(2.5) * smoothstep(px(4.0), px(2.5), dc) + srgb(50, 235, 240) * exp(-dc / px(8.0));
    return col;
}

vec3 shotBlur(vec2 p, float t) {
    vec2 c = ref(170.0, 290.0);
    vec2 q = p - c;
    float a = dot(q, JET), x = dot(q, vec2(-JET.y, JET.x));
    vec3 col = srgb(4, 5, 8);
    col += waveRings(p, ref(183.0, 316.0), JET, px(230.0), px(110.0), 2, srgb(18, 30, 48), px(18.0));
    col += srgb(140, 62, 32) * 0.75 * exp(-sq(a / px(70.0)) - sq(x / px(26.0))) * (0.8 + 0.4 * fbm(vec3(a * 12.0, x * 4.0, 1.0)));
    return col;
}

// --- Shot 3: gravitational waves ripple out as the jets drill out (7.15 - 8.07 s) ---------------

vec3 shotWaves(vec2 p, float t) {
    float s = t - 7.33;
    vec2 c = ref(180.0, 316.0);
    vec2 across = vec2(-JET.y, JET.x);
    float reach = px(117.0 + 100.0 * s);  // the heads fly apart
    vec2 h1 = c + JET * reach, h2 = c - JET * reach;
    float d1 = length(p - h1), d2 = length(p - h2);
    // A dark violet sky, lit magenta around the heads.
    vec3 col = srgb(13, 8, 20) + srgb(10, 6, 14) * band(-0.4, 0.4, -p.x + p.y);
    col += srgb(140, 30, 135) * (exp(-d1 / px(50.0)) + exp(-d2 / px(50.0))) * 0.22;
    col += starfield(p, 0.9);
    // Two rings, spreading: soft, lavender-blue.
    col += waveRings(p, c, JET, px(143.0 + 100.0 * s), px(135.0), 2, srgb(105, 120, 230) * 0.5, px(7.0));
    // The jet between the heads: an orange thread with knots, fuzzy.
    float along;
    float dl = segment(p, h2, h1, along);
    float knots = 0.65 + 0.35 * sin(along * 55.0 - t * 12.0);
    col += (srgb(255, 165, 75) * 1.1 * gauss(dl / px(4.0)) + srgb(235, 95, 45) * 0.6 * exp(-dl / px(11.0))) * knots;
    // The heads: magenta-white, clumpy, haloed.
    for (int i = 0; i < 2; ++i) {
        vec2 h = i == 0 ? h1 : h2;
        float d = length(p - h);
        vec2 q = (p - h) / px(38.0);
        float fluff = fbm(vec3(q * 2.4, t * 0.8 + float(i) * 7.0));
        float edge = px(36.0) * (0.75 + 0.5 * (fluff - 0.5));
        float body = smoothstep(edge, edge - px(14.0), d);
        float lumps = 0.6 + 0.8 * fbm(vec3(q * 5.0, 2.0 + float(i)));
        col += mix(srgb(235, 85, 225), srgb(255, 230, 255), saturate(lumps - 0.35)) * body * lumps * 1.05;
        col += srgb(235, 60, 215) * exp(-d / px(30.0)) * 0.9 + srgb(150, 30, 150) * exp(-d / px(70.0)) * 0.12;
    }
    // At the centre the remnant: a white spark between two fans of glowing gas.
    float dc = length(p - c);
    for (int i = -1; i <= 1; i += 2) {
        vec2 q = p - c - JET * px(19.0) * float(i);
        float e = length(vec2(dot(q, JET) / px(13.0), dot(q, across) / px(19.0)));
        float lobe = smoothstep(1.05, 0.55, e) * (0.45 + 0.55 * smoothstep(0.3, 0.95, e));
        col += srgb(185, 140, 255) * 0.85 * lobe;
    }
    col += vec3(2.5) * smoothstep(px(5.0), px(3.5), dc) + srgb(255, 230, 255) * exp(-dc / px(10.0)) * 0.8;
    return col;
}

// --- Shots 4-7: the burst, and the debris billowing into one cloud (8.07 - 10.9 s) --------------

// Orange smoke blown out both ways along the jet axis from c, widening with distance.
vec3 smokeJets(vec2 p, vec2 c, float t, float strength, float start) {
    vec2 across = vec2(-JET.y, JET.x);
    vec3 col = vec3(0.0);
    for (int i = -1; i <= 1; i += 2) {
        float a = dot(p - c, JET) * float(i);
        float x = dot(p - c, across);
        if (a <= 0.0) continue;
        float w = px(11.0) + 0.19 * a;
        float body = gauss(x / w) * band(start, start + px(30.0), a);
        float n1 = fbm(vec3(a * 9.0 - t * 0.9, x * 9.0, float(i) * 3.1));
        float n2 = fbm(vec3(a * 26.0 - t * 2.0, x * 26.0, float(i) * 5.7 + 1.0));
        float density = body * saturate(0.15 + 1.4 * n1 * n2 + 0.3 * n1);
        col += mix(srgb(130, 50, 15), srgb(225, 118, 38), saturate(n2 * 1.4 - 0.2)) * density;
    }
    return col * strength;
}

// The kilonova cloud's looks: cyan-white, then lavender, then grey smoke.
CloudLook blueCloud() {
    CloudLook l;
    l.deep = srgb(40, 110, 215);
    l.mid = srgb(145, 212, 255) * 1.15;
    l.lit = srgb(240, 252, 255) * 1.1;
    l.rim = srgb(205, 240, 255) * 1.3;
    l.grain = 0.45;
    return l;
}

CloudLook mixLook(CloudLook a, CloudLook b, float f) {
    CloudLook l;
    l.deep = mix(a.deep, b.deep, f);
    l.mid = mix(a.mid, b.mid, f);
    l.lit = mix(a.lit, b.lit, f);
    l.rim = mix(a.rim, b.rim, f);
    l.grain = mix(a.grain, b.grain, f);
    return l;
}

// The cocoons: hot gas the jet inflates on each side, soft, violet, brighter at a ragged rim,
// the jet's thread faintly pink through them.
vec3 cocoons(vec2 p, vec2 c, float offset, float size, float fade) {
    vec2 across = vec2(-JET.y, JET.x);
    vec3 col = vec3(0.0);
    for (int i = -1; i <= 1; i += 2) {
        vec2 q = p - c - JET * offset * float(i);
        float a = dot(q, JET), x = dot(q, across);
        float wob = fbm(vec3(q / size * 2.0, float(i) * 4.0));
        float e = length(vec2(a / (size * 0.85), x / (size * 1.05))) * (0.88 + 0.3 * (wob - 0.5));
        float lobe = smoothstep(1.05, 0.65, e) * (0.35 + 0.65 * smoothstep(0.4, 0.97, e));
        float streaks = 0.6 + 0.7 * fbm(vec3(a / size * 3.0, x / size * 1.2, float(i) + 9.0));
        col += srgb(100, 78, 255) * 1.1 * lobe * streaks;
        col += srgb(235, 90, 200) * gauss(x / (size * 0.22)) * smoothstep(1.05, 0.5, e) * 0.45;
    }
    return col * fade;
}

vec3 shotBurst(vec2 p, float t) {
    float s = saturate((t - 8.13) / 0.8);
    vec2 c = ref(179.0, 322.0);
    vec3 col = srgb(7, 4, 3) + starfield(p, 0.8);
    // The rings and the heads' glow linger for a moment from the last shot.
    float linger = 1.0 - band(T_BURST, T_BURST + 0.3, t);
    col += waveRings(p, c, JET, px(260.0), px(135.0), 2, srgb(22, 28, 66), px(10.0)) * linger;
    col += srgb(200, 45, 180) * (exp(-length(p - ref(0.0, 215.0)) / px(45.0)) + exp(-length(p - ref(385.0, 440.0)) / px(45.0))) *
           linger * 0.45;
    col += smokeJets(p, c, t, 1.0, px(18.0));
    float remnant = px(20.0 + 25.0 * s);
    float size = px(16.0 + 32.0 * s);
    col += cocoons(p, c, remnant * 0.85 + size * 0.6, size, 1.0);
    // The remnant: a flattened grey-white disk of debris that puffs up into a ball.
    Cloud cl = cloudAt(c, remnant);
    cl.axis = normalize(vec3(JET, 0.35));
    cl.squash = mix(0.55, 0.92, s);
    cl.seed = 3.0;
    cl.spin = 0.3 * t;
    cl.rough = mix(0.1, 0.28, s);
    cl.scale = 4.5;
    CloudLook grey;
    grey.deep = srgb(110, 118, 140);
    grey.mid = srgb(190, 196, 210);
    grey.lit = srgb(240, 242, 248);
    grey.rim = srgb(180, 190, 225) * 0.7;
    grey.grain = 0.4;
    vec4 cloud = drawCloud(viewRay(p), cl, mixLook(grey, blueCloud(), s * s));
    return mix(col, cloud.rgb, cloud.a);
}

// The cloud's centre and radius (frame) from 9 s to the end of the "heavy elements" shot,
// matched to its outline in the reference.
void kilonovaCloud(float t, out vec2 c, out float r) {
    if (t < T_WHIP) {
        float u = (t - 9.07) / 0.8;
        c = mix(ref(181.0, 319.0), ref(174.0, 311.0), saturate(u));
        r = px(50.0 + 34.0 * u);
    } else {
        c = mix(ref(170.0, 312.0), ref(166.0, 305.0), saturate((t - 10.13) / 0.77));
        r = px(128.0 + 55.0 * (t - 10.4));
    }
}

// The shots around the cloud: its background (sky, halo, jets behind it) and, separately, the
// cloud itself (colour and coverage), for the caller to lay over it.
vec3 cloudShot(vec2 p, float t, float jetStrength, CloudLook look, float halo, out vec4 cloud) {
    vec2 c;
    float r;
    kilonovaCloud(t, c, r);
    float d = length(p - c);
    vec3 col = srgb(8, 5, 5) + starfield(p, 1.0);
    col += srgb(12, 15, 38) * exp(-d / 0.3) * halo;
    col += smokeJets(p, c, t, jetStrength, px(10.0));
    col += srgb(120, 175, 255) * exp(-max(d - r, 0.0) / px(15.0)) * 0.6 * halo;
    Cloud cl = cloudAt(c, r);
    cl.seed = 3.0;
    cl.spin = 0.3 * t;
    cl.squash = 0.92;
    cl.axis = normalize(vec3(JET, 0.35));
    cloud = drawCloud(viewRay(p), cl, look);
    return col;
}

vec3 shotMerge(vec2 p, float t) {
    vec4 cloud;
    vec3 col = cloudShot(p, t, 0.7, blueCloud(), 1.0, cloud);
    // The cocoons' far rims, flying off as fuzzy violet arcs (the halves facing away).
    vec2 c = ref(181.0, 319.0);
    float fade = 1.0 - band(9.3, 9.8, t);
    for (int i = -1; i <= 1; i += 2) {
        vec2 a = c + JET * px(55.0 + 60.0 * (t - 9.07)) * float(i);
        float along = dot(p - a, JET) * float(i);
        float wob = (fbm(vec3((p - a) * 40.0, float(i))) - 0.5) * px(8.0);
        float d = ellipseLine(p, a, JET, px(46.0), px(62.0)) + wob;
        float line = gauss(d / px(4.0)) + 0.4 * gauss(d / px(12.0)) + 0.12 * step(d, 0.0);
        col += srgb(120, 85, 255) * line * band(-px(10.0), px(25.0), along) * fade * 1.2;
    }
    return mix(col, cloud.rgb, cloud.a);
}

vec3 shotWhip(vec2 p, float t) {
    vec2 c = ref(150.0, 268.0);
    float r = px(80.0);
    vec2 q = p - c;
    float shape = exp(-sq(sq(q.x / r))) * exp(-sq(sq(q.y / (r * 1.1))));
    float streaks = 0.55 + 0.9 * fbm(vec3(q.x * 90.0, 1.0, 2.0));
    return srgb(6, 6, 10) + srgb(140, 190, 245) * shape * streaks * 1.1 + srgb(25, 30, 70) * exp(-length(q) / 0.25) * 0.5;
}

vec3 shotHeavy(vec2 p, float t) {
    vec4 cloud;
    vec3 col = cloudShot(p, t, 0.45, blueCloud(), 1.0, cloud);
    return mix(col, cloud.rgb, cloud.a);
}

// --- Shots 8-11: the kilonova's glow, and a new black hole (10.9 - 14.72 s) ---------------------

// The red remnant at the upper left: a clumpy ball, pink at first, darkening to red, drifting
// off the frame's edge.
vec2 ballCentre(float t) { return mix(ref(20.0, 160.0), ref(-4.0, 140.0), band(11.3, 14.6, t)); }

// The shell's colours: orange, then (11.95 - 12.9 s) olive, teal and blue; tl is the local
// time, later toward the rim. deep: its body, bright: the rim and knot.
void shellColours(float tl, out vec3 deep, out vec3 bright) {
    vec3 d0 = srgb(170, 82, 6), b0 = srgb(255, 185, 35);
    vec3 d1 = srgb(115, 108, 45), b1 = srgb(205, 195, 80);
    vec3 d2 = srgb(35, 105, 150), b2 = srgb(95, 225, 215);
    vec3 d3 = srgb(20, 115, 220), b3 = srgb(60, 215, 255);
    vec3 d4 = srgb(25, 135, 245), b4 = srgb(95, 250, 255);
    float k = tl;
    deep = k < 12.15 ? mix(d0, d1, band(11.92, 12.15, k))
         : k < 12.45 ? mix(d1, d2, band(12.15, 12.45, k))
         : k < 12.8  ? mix(d2, d3, band(12.45, 12.8, k))
                     : mix(d3, d4, band(13.0, 14.3, k));
    bright = k < 12.15 ? mix(b0, b1, band(11.92, 12.15, k))
           : k < 12.45 ? mix(b1, b2, band(12.15, 12.45, k))
           : k < 12.8  ? mix(b2, b3, band(12.45, 12.8, k))
                       : mix(b3, b4, band(13.0, 14.3, k));
}

// The glowing shell of ejecta: a sphere of fibrous gas seen from inside its cone of light,
// only its far half (lower right, away from the remnant) shining, its fibres streaming
// outward, its rim there sharp and bright, a bright knot where the remnant's haze meets it.
// Returns the light and how much of what lies behind it (the haze) it hides.
vec4 shell(vec2 p, float t) {
    float grow = band(11.2, 14.7, t);
    vec2 c = mix(ref(222.0, 470.0), ref(232.0, 462.0), grow);
    float R = px(mix(158.0, 222.0, grow));
    vec2 q = (p - c) / R;
    float rho = length(q);
    if (rho > 1.15) return vec4(0.0);
    vec2 knot = mix(ref(238.0, 452.0), ref(240.0, 385.0), band(12.3, 13.2, t));
    float dk = length(p - knot);
    vec3 deep, bright;
    shellColours(t - 0.25 * rho, deep, bright);
    float z = sqrt(max(1.0 - rho * rho, 0.0));
    float limb = 1.0 / (z + 0.2);
    // Fibres: noise stretched along the radius (fine around, coarse outward), mottled.
    float ang = atan(q.y, q.x);
    float fibres = fbm(vec3(cos(ang) * 9.0, sin(ang) * 9.0, rho * 2.5 - t * 0.04)) *
                   fbm(vec3(cos(ang) * 22.0 + 3.0, sin(ang) * 22.0, rho * 5.0));
    float mottle = fbm(vec3(q * 7.0, 1.0 + t * 0.05));
    float fil = 0.15 + 2.2 * fibres * (0.6 + 0.8 * mottle);
    // The far side shines; toward the remnant it thins away.
    vec2 lr = normalize(vec2(0.55, -0.83));
    float side = dot(q, lr);
    float lit = band(-0.75, 0.55, side);
    float edge = smoothstep(1.07, 0.96, rho);
    vec3 col = deep * fil * limb * 0.34 * lit * edge;
    // The bright rim on the far side, fuzzy with fibres.
    float dir = dot(q, lr) / max(rho, 1e-3);
    float rim = gauss((rho - 0.985 - 0.02 * (fibres - 0.3)) / 0.022) * band(-0.05, 0.6, dir);
    col += bright * rim * (0.5 + 1.2 * fibres) * 1.1;
    // The glowing knot.
    col += bright * exp(-dk / px(65.0)) * 0.9 * edge * (0.6 + 0.6 * mottle);
    // A faint glow just outside the rim.
    col += deep * exp(-max(rho - 1.0, 0.0) / 0.035) * step(1.0, rho) * 0.25 * band(0.0, 0.6, dir);
    float opaque = mix(0.3, 0.95, band(12.0, 12.6, t));
    return vec4(col, edge * opaque);
}

vec3 remnantBall(vec2 p, float t, vec3 col, float show) {
    vec2 c = ballCentre(t);
    float r = px(mix(25.0, 31.0, band(11.3, 14.6, t)));
    Cloud cl = cloudAt(c, r);
    cl.seed = 41.0;
    cl.spin = 0.2 * t;
    cl.rough = 0.22;
    cl.scale = 5.5;
    CloudLook pink, red;
    pink.deep = srgb(150, 60, 85);
    pink.mid = srgb(225, 140, 155);
    pink.lit = srgb(255, 225, 230);
    pink.rim = srgb(240, 170, 180) * 0.5;
    pink.grain = 0.6;
    red.deep = srgb(70, 14, 6);
    red.mid = srgb(160, 45, 18);
    red.lit = srgb(215, 95, 55);
    red.rim = srgb(170, 55, 25) * 0.4;
    red.grain = 0.6;
    float darken = band(11.75, 12.1, t);
    vec4 b = drawCloud(viewRay(p), cl, mixLook(pink, red, darken));
    float d = length(p - c);
    col += mix(srgb(200, 90, 100), srgb(160, 40, 10), darken) * exp(-max(d - r, 0.0) / px(16.0)) * 0.4 * show;
    return mix(col, b.rgb, b.a * show);
}

vec3 shotShell(vec2 p, float t, float shellShow, float ballShow) {
    vec3 col = srgb(4, 4, 5) + starfield(p, 1.0);
    // The orange haze: a cone of light from the remnant spreading to the shell.
    vec2 ball = ballCentre(t);
    vec2 target = ref(225.0, 460.0);
    float along;
    float d = segment(p, ball, target, along);
    float w = mix(px(10.0), px(130.0), along);
    float haze = gauss(d / w) * (0.2 + 0.8 * along) * (0.75 + 0.5 * fbm(vec3(p * 10.0, 2.0)));
    vec3 hazeCol = srgb(185, 85, 22) * haze * 0.75;
    // Later the haze gathers into a glow at the shell's upper-left edge.
    hazeCol += srgb(170, 65, 10) * exp(-length(p - ref(110.0, 300.0)) / px(65.0)) * 0.6 * band(12.4, 13.2, t);
    vec4 sh = shell(p, t);
    col += hazeCol * (1.0 - sh.a * shellShow);
    col += sh.rgb * shellShow;
    col = remnantBall(p, t, col, ballShow);
    col += starfield(p, 0.5);  // the reference's stars show over the glow
    return col;
}

vec3 shotFade(vec2 p, float t) {
    // The cloud turns lavender, then grey, and thins to smoke over the glowing shell.
    float lav = band(T_FADE, 11.03, t);
    float smoke = band(11.08, 11.25, t);
    CloudLook lavender;
    lavender.deep = srgb(110, 105, 175);
    lavender.mid = srgb(185, 182, 235);
    lavender.lit = srgb(240, 238, 252);
    lavender.rim = srgb(210, 208, 250) * 1.0;
    lavender.grain = 0.4;
    CloudLook grey;
    grey.deep = srgb(50, 48, 56);
    grey.mid = srgb(100, 97, 106);
    grey.lit = srgb(140, 136, 142);
    grey.rim = srgb(110, 105, 120) * 0.4;
    grey.grain = 0.7;
    CloudLook look = mixLook(mixLook(blueCloud(), lavender, lav), grey, smoke);
    vec4 cloud;
    vec3 bg = cloudShot(p, t, 0.3 * (1.0 - smoke), look, 1.0 - smoke, cloud);
    // The shell's scene comes up behind as the cloud's own background goes.
    float under = band(11.05, 11.2, t);
    bg = mix(bg, shotShell(p, t, band(11.1, 11.3, t), band(10.98, 11.15, t)), under);
    float alpha = cloud.a * saturate(1.0 - smoke * (0.4 + 0.9 * fbm(vec3(p * 9.0, t))));
    return mix(bg, cloud.rgb, alpha);
}

// --- The last shot: the two stars again, calm, before it all (14.9 - 17.6 s) ---------------------

vec3 shotCalm(vec2 p, float t) {
    vec2 c = ref(172.0, 322.0);
    vec3 bg = srgb(2, 4, 8) + srgb(8, 40, 62) * exp(-dot(p - c, p - c) / 0.06) + starfield(p, 1.0);
    vec3 col = binaryStars(p, bg, 0.75);
    // Fading in from black.
    return col * pow(band(T_CALM + 0.05, 15.6, t), 2.0);
}

// --- Main --------------------------------------------------------------------------------------------

void main() {
    pixel = 1.0 / uFrameHeight;
    vec2 p = (gl_FragCoord.xy - 0.5 * uResolution) / uFrameHeight;
    float t = uTime;
    vec3 col;
    if (t < T_FLASH) col = shotInspiral(p, t);
    else if (t < T_WARP) col = shotFlash(p, t);
    else if (t < T_STREAK) col = shotWarp(p, t);
    else if (t < T_BLUR) col = shotStreak(p, t);
    else if (t < T_WAVES) col = shotBlur(p, t);
    else if (t < T_BURST) col = shotWaves(p, t);
    else if (t < T_MERGE) col = shotBurst(p, t);
    else if (t < T_WHIP) col = shotMerge(p, t);
    else if (t < T_HEAVY) col = shotWhip(p, t);
    else if (t < T_FADE) col = shotHeavy(p, t);
    else if (t < T_SHELL) col = shotFade(p, t);
    else if (t < T_BLACK) col = shotShell(p, t, 1.0, 1.0) * (1.0 - band(14.6, T_BLACK, t));
    else if (t < T_CALM) col = vec3(0.0);
    else col = shotCalm(p, t);
    outColor = vec4(max(col, 0.0), 1.0);
}
