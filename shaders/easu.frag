#version 460 core

// Edge-adaptive spatial upsampling, after AMD's FidelityFX Super Resolution 1 (EASU, MIT
// licence): the tone-mapped image, traced below the window's resolution, is enlarged with
// a Lanczos-2-like kernel stretched along the local edge and squeezed across it, so edges
// stay crisp instead of blurring as bilinear filtering blurs them.
//
// For each output pixel, 12 input pixels around it (a 4x4 block without its corners):
//          b c
//        e f g h
//        i j k l
//          n o
// The direction and strength of the edge come from luma differences around the four
// nearest (f, g, j, k), weighted bilinearly. The kernel, a polynomial approximation of
// Lanczos 2, is rotated to the edge, stretched along it (up to sqrt 2) and narrowed across
// it, its negative lobe weaker where there is no edge. The result is clamped to the four
// nearest pixels' range, so it cannot ring.

in vec2 vUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D uInput;
uniform vec2 uInputSize;
uniform vec2 uOutputSize;

vec3 tapAt(ivec2 p) {
    return texelFetch(uInput, clamp(p, ivec2(0), ivec2(uInputSize) - 1), 0).rgb;
}

// FSR's luma, twice the usual weights' scale: 0.5 R + G + 0.5 B.
float lumaOf(vec3 c) {
    return 0.5 * c.r + c.g + 0.5 * c.b;
}

// Adds one of the four nearest pixels' edge direction and length, with bilinear weight w,
// from its neighbours: A above, B left, C itself, D right, E below.
void edgeOf(inout vec2 dir, inout float len, float w, float lA, float lB, float lC, float lD, float lE) {
    float lenX = max(abs(lD - lC), abs(lC - lB));
    lenX = lenX > 0.0 ? 1.0 / lenX : 0.0;
    float dirX = lD - lB;
    dir.x += dirX * w;
    lenX = clamp(abs(dirX) * lenX, 0.0, 1.0);
    len += lenX * lenX * w;

    float lenY = max(abs(lE - lC), abs(lC - lA));
    lenY = lenY > 0.0 ? 1.0 / lenY : 0.0;
    float dirY = lE - lA;
    dir.y += dirY * w;
    lenY = clamp(abs(dirY) * lenY, 0.0, 1.0);
    len += lenY * lenY * w;
}

// One tap of the rotated, stretched kernel at offset off from the sample point.
void addTap(inout vec3 sum, inout float weight, vec2 off, vec2 dir, vec2 len2, float lob, float clp, vec3 c) {
    vec2 v = vec2(off.x * dir.x + off.y * dir.y, -off.x * dir.y + off.y * dir.x) * len2;
    float d2 = min(dot(v, v), clp);
    // (25/16 (2/5 d2 - 1)^2 - (25/16 - 1)) (lob d2 - 1)^2: a window times a lobe.
    float wB = 0.4 * d2 - 1.0;
    float wA = lob * d2 - 1.0;
    wB *= wB;
    wA *= wA;
    wB = 1.5625 * wB - 0.5625;
    float w = wB * wA;
    sum += c * w;
    weight += w;
}

void main() {
    // Where this output pixel's center falls among the input pixels.
    vec2 pp = gl_FragCoord.xy * uInputSize / uOutputSize - 0.5;
    vec2 fp = floor(pp);
    pp -= fp;
    ivec2 p = ivec2(fp);

    vec3 b = tapAt(p + ivec2(0, -1)), c = tapAt(p + ivec2(1, -1));
    vec3 e = tapAt(p + ivec2(-1, 0)), f = tapAt(p), g = tapAt(p + ivec2(1, 0)), h = tapAt(p + ivec2(2, 0));
    vec3 i = tapAt(p + ivec2(-1, 1)), j = tapAt(p + ivec2(0, 1)), k = tapAt(p + ivec2(1, 1)), l = tapAt(p + ivec2(2, 1));
    vec3 n = tapAt(p + ivec2(0, 2)), o = tapAt(p + ivec2(1, 2));
    float bL = lumaOf(b), cL = lumaOf(c), eL = lumaOf(e), fL = lumaOf(f), gL = lumaOf(g), hL = lumaOf(h);
    float iL = lumaOf(i), jL = lumaOf(j), kL = lumaOf(k), lL = lumaOf(l), nL = lumaOf(n), oL = lumaOf(o);

    vec2 dir = vec2(0.0);
    float len = 0.0;
    edgeOf(dir, len, (1.0 - pp.x) * (1.0 - pp.y), bL, eL, fL, gL, jL);
    edgeOf(dir, len, pp.x * (1.0 - pp.y), cL, fL, gL, hL, kL);
    edgeOf(dir, len, (1.0 - pp.x) * pp.y, fL, iL, jL, kL, nL);
    edgeOf(dir, len, pp.x * pp.y, gL, jL, kL, lL, oL);

    // Normalize the direction (no edge: along x); the length from 0 (none) to 1 (strong).
    float dirR = dot(dir, dir);
    bool noEdge = dirR < 1.0 / 32768.0;
    dirR = noEdge ? 1.0 : inversesqrt(dirR);
    dir.x = noEdge ? 1.0 : dir.x;
    dir *= dirR;
    len *= 0.5;
    len *= len;
    // Stretch along the edge (1 on an axis, sqrt 2 on a diagonal), narrow across it, and
    // soften the negative lobe where the edge is weak.
    float stretch = dot(dir, dir) / max(abs(dir.x), abs(dir.y));
    vec2 len2 = vec2(1.0 + (stretch - 1.0) * len, 1.0 - 0.5 * len);
    float lob = 0.5 + ((1.0 / 4.0 - 0.04) - 0.5) * len;
    float clp = 1.0 / lob;

    vec3 sum = vec3(0.0);
    float weight = 0.0;
    addTap(sum, weight, vec2(0.0, -1.0) - pp, dir, len2, lob, clp, b);
    addTap(sum, weight, vec2(1.0, -1.0) - pp, dir, len2, lob, clp, c);
    addTap(sum, weight, vec2(-1.0, 0.0) - pp, dir, len2, lob, clp, e);
    addTap(sum, weight, vec2(0.0, 0.0) - pp, dir, len2, lob, clp, f);
    addTap(sum, weight, vec2(1.0, 0.0) - pp, dir, len2, lob, clp, g);
    addTap(sum, weight, vec2(2.0, 0.0) - pp, dir, len2, lob, clp, h);
    addTap(sum, weight, vec2(-1.0, 1.0) - pp, dir, len2, lob, clp, i);
    addTap(sum, weight, vec2(0.0, 1.0) - pp, dir, len2, lob, clp, j);
    addTap(sum, weight, vec2(1.0, 1.0) - pp, dir, len2, lob, clp, k);
    addTap(sum, weight, vec2(2.0, 1.0) - pp, dir, len2, lob, clp, l);
    addTap(sum, weight, vec2(0.0, 2.0) - pp, dir, len2, lob, clp, n);
    addTap(sum, weight, vec2(1.0, 2.0) - pp, dir, len2, lob, clp, o);

    // No ringing: within the four nearest pixels' range.
    vec3 lo = min(min(f, g), min(j, k)), hi = max(max(f, g), max(j, k));
    outColor = vec4(clamp(sum / max(weight, 1e-6), lo, hi), 1.0);
}
