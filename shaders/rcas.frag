#version 460 core

// Robust contrast-adaptive sharpening, after AMD's FidelityFX Super Resolution 1 (RCAS,
// MIT licence): after the upscale, each pixel is sharpened against its four neighbours by
// a negative lobe as strong as it can be without pushing the result outside the
// neighbourhood's range (so it cannot clip or ring), at most -(1/4 - 1/16), scaled by the
// sharpness. Then the dither that hides banding in the dark sky.

in vec2 vUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D uInput;  // the upscaled image, at the output size
uniform float uSharpness;  // 2^-stops: 1 is the sharpest

const float kLimit = 0.25 - 1.0 / 16.0;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 last = textureSize(uInput, 0) - 1;
    //   b
    // d e f
    //   h
    vec3 b = texelFetch(uInput, clamp(p + ivec2(0, 1), ivec2(0), last), 0).rgb;
    vec3 d = texelFetch(uInput, clamp(p + ivec2(-1, 0), ivec2(0), last), 0).rgb;
    vec3 e = texelFetch(uInput, p, 0).rgb;
    vec3 f = texelFetch(uInput, clamp(p + ivec2(1, 0), ivec2(0), last), 0).rgb;
    vec3 h = texelFetch(uInput, clamp(p + ivec2(0, -1), ivec2(0), last), 0).rgb;
    vec3 mn = min(min(min(b, d), min(f, h)), e);
    vec3 mx = max(max(max(b, d), max(f, h)), e);
    // The strongest lobe that keeps the result above the least and below the most.
    vec3 hitMin = mn / (4.0 * mx + 1e-6);
    vec3 hitMax = (1.0 - mx) / (4.0 * mn - 4.0 - 1e-6);
    vec3 lobes = max(-hitMin, hitMax);
    float lobe = max(-kLimit, min(max(lobes.r, max(lobes.g, lobes.b)), 0.0)) * uSharpness;
    vec3 c = (lobe * (b + d + f + h) + e) / (4.0 * lobe + 1.0);
    // Interleaved gradient noise dither hides banding in the dark sky.
    float noise = fract(52.9829189 * fract(dot(gl_FragCoord.xy, vec2(0.06711056, 0.00583715))));
    outColor = vec4(c + (noise - 0.5) / 255.0, 1.0);
}
