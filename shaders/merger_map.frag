#version 460 core

// The neutron star merger, simulation view: the matter seen from above, as numerical
// relativity papers show it. Each pixel is the column density (g/cm^2) through the
// orbital plane, on a logarithmic scale through the "inferno" color map (black, purple,
// orange, yellow), the hottest matter a little whiter. Written straight to the 8-bit
// image (no tone mapping), with a faint grid every uGridStep.
//
// Once the hydrodynamics has stopped, the view turns edge-on and shows what follows as a
// slice through the spin axis: the kilonova's ejecta (src/merger.h; the same shapes the
// cinematic view draws), bright where dense and colored by what they are made of, the
// lanthanide-poor (blue) ejecta toward the poles and the lanthanide-rich (red) ones in the
// orbital plane, with the funnel the jet clears and its cocoon.

layout(location = 0) out vec4 outColor;

#include "noise.glsl"

layout(binding = 0) uniform sampler2D uMap;  // log10 column density, log10 temperature
uniform int uEdgeOn;          // the slice through the axis instead of the map
uniform vec4 uEjecta[3];      // per component: outer radius (simulation units), log10 mean density (g/cm^3)
uniform vec4 uFunnel;         // the jet: head distance, half-opening angle, cocoon radius, breakout distance

uniform vec2 uViewport;     // lower-left pixel of this view
uniform vec2 uResolution;   // its size
uniform vec2 uCenter;       // what the view is centred on (simulation units)
uniform float uHalfHeight;  // half its height (simulation units)
uniform vec2 uMapMin;
uniform float uMapSize;
uniform vec2 uLogRange;     // column densities mapped to the ends of the color map
uniform float uGridStep;    // simulation units between grid lines

// Matplotlib's "inferno" (van der Walt & Smith), polynomial fit by Matt Zucker (CC0).
vec3 inferno(float t) {
    const vec3 c0 = vec3(0.0002189403691192265, 0.001651004631001012, -0.01948089843709184);
    const vec3 c1 = vec3(0.1065134194856116, 0.5639564367884091, 3.932712388889277);
    const vec3 c2 = vec3(11.60249308247187, -3.972853965665698, -15.9423941062914);
    const vec3 c3 = vec3(-41.70399613139459, 17.43639888205313, 44.35414519872813);
    const vec3 c4 = vec3(77.162935699427, -33.40235894210092, -81.80730925738993);
    const vec3 c5 = vec3(-71.31942824499214, 32.62606426397723, 73.20951985803202);
    const vec3 c6 = vec3(25.13112622477341, -12.24266895238567, -23.07032500287172);
    t = clamp(t, 0.0, 1.0);
    return clamp(c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6))))), 0.0, 1.0);
}

// The ejecta at (x, 0, z): log10 of the density (g/cm^3) and the color of what is there,
// mixed by mass from the three components (blue, purple, red), with the same shapes as the
// cinematic view (merger_view.frag) without its finest noise.
vec3 componentColor(int k) {
    return k == 0 ? vec3(0.25, 0.55, 1.0) : (k == 1 ? vec3(0.85, 0.35, 0.75) : vec3(1.0, 0.45, 0.08));
}

float ejecta(vec2 p, out vec3 hue) {
    float r = length(p);
    float c = abs(p.y) / max(r, 1e-6);
    float rmax = max(uEjecta[0].x * 1.15, max(uEjecta[1].x * 1.15, uEjecta[2].x * 1.6));
    vec3 q = vec3(p.x, 0.0, p.y) / max(rmax, 1e-6);
    float lumps = 0.5 + 0.9 * fbmFiltered(q * 4.0, 5, 0.0) + 0.5 * ridgedFiltered(q * 9.0 + 7.3, 4, 0.0);
    float rho = 0.0;
    hue = vec3(0.0);
    for (int k = 0; k < 3; ++k) {
        if (uEjecta[k].x <= 0.0) continue;
        vec2 ab = k == 0 ? vec2(0.8, 1.15) : (k == 1 ? vec2(1.15, 0.85) : vec2(1.6, 0.55));
        float stretch = inversesqrt((1.0 - c * c) / (ab.x * ab.x) + c * c / (ab.y * ab.y));
        float s = r / (uEjecta[k].x * stretch);
        if (s >= 1.0) continue;
        float ang = k == 0 ? 0.06 + 1.8 * smoothstep(0.3, 0.85, c)
                  : (k == 1 ? 0.45 + 0.55 * (1.0 - c * c) : 0.15 + 2.4 * exp(-c * c / 0.09));
        float shape = (1.0 - 0.6 * s * s) * smoothstep(1.0, 0.85, s) * ang * smoothstep(0.15, 0.85, lumps);
        float rk = pow(10.0, uEjecta[k].y) * shape;
        rho += rk;
        hue += rk * componentColor(k);
    }
    hue /= max(rho, 1e-40);
    // The jet's funnel, and the cocoon's shell around where it broke out.
    if (uFunnel.x > 0.0) {
        float th = acos(c);
        if (th < uFunnel.y * 1.6 && r < uFunnel.x) rho *= 0.03;
        if (uFunnel.z > 0.0) {
            float d = length(vec2(p.x, abs(p.y) - uFunnel.w));
            float shell = exp(-pow((d - uFunnel.z) / (0.12 * uFunnel.z), 2.0));
            rho *= 1.0 + 4.0 * shell;
            hue = mix(hue, vec3(0.75, 0.5, 1.0), shell);
        }
    }
    return log(max(rho, 1e-40)) / log(10.0);
}

void main() {
    vec2 uv = (gl_FragCoord.xy - uViewport) / uResolution;
    float aspect = uResolution.x / uResolution.y;
    vec2 p = uCenter + (uv * 2.0 - 1.0) * vec2(aspect, 1.0) * uHalfHeight;
    vec2 m = (p - uMapMin) / uMapSize;
    vec3 c = vec3(0.0);
    if (uEdgeOn == 1) {
        // Brightness from the density (log scale), color from the composition.
        vec3 hue;
        float t = (ejecta(p, hue) - uLogRange.x) / (uLogRange.y - uLogRange.x);
        if (t > 0.0) c = hue * (0.15 + 1.1 * pow(clamp(t, 0.0, 1.0), 1.5)) + vec3(0.35) * pow(clamp(t, 0.0, 1.0), 6.0);
    } else if (all(greaterThan(m, vec2(0.0))) && all(lessThan(m, vec2(1.0)))) {
        vec2 s = texture(uMap, m).rg;
        float t = (s.x - uLogRange.x) / (uLogRange.y - uLogRange.x);
        c = inferno(t);
        // Shocked, hot matter glows whiter.
        float heat = clamp((s.y - 10.0) / 1.2, 0.0, 1.0) * smoothstep(0.55, 0.9, t);
        c = mix(c, vec3(1.0, 0.97, 0.85), 0.35 * heat);
    }
    // A faint grid, as on the plots of simulations.
    vec2 g = p / uGridStep;
    vec2 w = fwidth(g);
    vec2 line = 1.0 - smoothstep(vec2(0.0), w * 1.2, abs(fract(g + 0.5) - 0.5));
    c = max(c, vec3(0.13, 0.13, 0.16) * max(line.x, line.y));
    // sRGB-ish encoding: the color map is defined in display values already.
    outColor = vec4(c, 1.0);
}
