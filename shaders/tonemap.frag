#version 460 core

// Combines the HDR image with bloom, applies exposure and a tone curve, and
// encodes to sRGB.

in vec2 vUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D uHdr;
layout(binding = 1) uniform sampler2D uBloom;
uniform float uExposure;
uniform float uBloomStrength;
uniform float uBloomLevels;  // the bloom chain sums this many levels
uniform int uTonemapper;     // 0 = ACES, 1 = Reinhard, 2 = clamp

// Narkowicz's fit of the ACES filmic curve.
vec3 aces(vec3 x) {
    return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}

vec3 linearToSrgb(vec3 c) {
    return mix(c * 12.92, 1.055 * pow(c, vec3(1.0 / 2.4)) - 0.055, step(0.0031308, c));
}

void main() {
    vec3 hdr = texture(uHdr, vUV).rgb;
    vec3 bloom = texture(uBloom, vUV).rgb / uBloomLevels;
    vec3 c = mix(hdr, bloom, uBloomStrength) * uExposure;

    if (uTonemapper == 0) c = aces(c);
    else if (uTonemapper == 1) c = c / (1.0 + c);
    else c = clamp(c, 0.0, 1.0);

    c = linearToSrgb(c);
    // Interleaved gradient noise dither hides banding in the dark sky.
    float noise = fract(52.9829189 * fract(dot(gl_FragCoord.xy, vec2(0.06711056, 0.00583715))));
    outColor = vec4(c + (noise - 0.5) / 255.0, 1.0);
}
