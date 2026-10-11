#version 460 core

// 13-tap downsample (Jimenez, "Next Generation Post Processing in Call of Duty:
// Advanced Warfare", 2014). The first pass uses a Karis average to keep single
// very bright pixels from turning into flickering blobs.

in vec2 vUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D uSrc;
uniform vec2 uSrcTexel;
uniform int uFirst;

vec3 tap(float x, float y) {
    return texture(uSrc, vUV + vec2(x, y) * uSrcTexel).rgb;
}

float luma(vec3 c) {
    return dot(c, vec3(0.2126, 0.7152, 0.0722));
}

void main() {
    vec3 a = tap(-2.0, 2.0), b = tap(0.0, 2.0), c = tap(2.0, 2.0);
    vec3 d = tap(-2.0, 0.0), e = tap(0.0, 0.0), f = tap(2.0, 0.0);
    vec3 g = tap(-2.0, -2.0), h = tap(0.0, -2.0), i = tap(2.0, -2.0);
    vec3 j = tap(-1.0, 1.0), k = tap(1.0, 1.0), l = tap(-1.0, -1.0), m = tap(1.0, -1.0);

    vec3 result;
    if (uFirst == 1) {
        vec3 g0 = (j + k + l + m) * 0.25;
        vec3 g1 = (a + b + d + e) * 0.25;
        vec3 g2 = (b + c + e + f) * 0.25;
        vec3 g3 = (d + e + g + h) * 0.25;
        vec3 g4 = (e + f + h + i) * 0.25;
        float w0 = 0.5 / (1.0 + luma(g0));
        float w1 = 0.125 / (1.0 + luma(g1));
        float w2 = 0.125 / (1.0 + luma(g2));
        float w3 = 0.125 / (1.0 + luma(g3));
        float w4 = 0.125 / (1.0 + luma(g4));
        result = (g0 * w0 + g1 * w1 + g2 * w2 + g3 * w3 + g4 * w4) / (w0 + w1 + w2 + w3 + w4);
    } else {
        result = e * 0.125 + (a + c + g + i) * 0.03125 + (b + d + f + h) * 0.0625 + (j + k + l + m) * 0.125;
    }
    outColor = vec4(max(result, vec3(0.0)), 1.0);
}
