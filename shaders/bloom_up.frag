#version 460 core

// 3x3 tent-filter upsample, blended additively onto the next larger level.

in vec2 vUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D uSrc;
uniform vec2 uSrcTexel;

vec3 tap(float x, float y) {
    return texture(uSrc, vUV + vec2(x, y) * uSrcTexel).rgb;
}

void main() {
    vec3 sum = tap(0.0, 0.0) * 4.0;
    sum += (tap(0.0, 1.0) + tap(-1.0, 0.0) + tap(1.0, 0.0) + tap(0.0, -1.0)) * 2.0;
    sum += tap(-1.0, 1.0) + tap(1.0, 1.0) + tap(-1.0, -1.0) + tap(1.0, -1.0);
    outColor = vec4(sum / 16.0, 1.0);
}
