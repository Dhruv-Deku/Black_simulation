#version 460 core

// Bakes the Sun's 30.4 nm chromosphere (shaders/sdo.glsl) into an equirectangular map of
// its surface, in the rotating frame: u east longitude, v latitude (south pole at 0).

in vec2 vUV;
layout(location = 0) out vec4 outColor;

uniform int uRegionCount;
uniform vec4 uRegions[8];   // active regions: direction, angular radius
uniform float uGranuleFreq; // convection cells per radian (the network's are 25 times larger)
uniform float uTime;        // the chromosphere's evolution
uniform vec2 uMapSize;

const float PI = 3.14159265358979;
#include "noise.glsl"
#include "sdo.glsl"

void main() {
    float lon = (vUV.x - 0.5) * 2.0 * PI;
    float lat = (vUV.y - 0.5) * PI;
    vec3 n = vec3(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat));
    // Filtered to a texel's size on the surface (its height: the width shrinks poleward).
    outColor = vec4(sdoIntensity(n, PI / uMapSize.y, uTime), 0.0, 0.0, 1.0);
}
