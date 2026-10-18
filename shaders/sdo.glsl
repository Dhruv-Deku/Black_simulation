// The Sun as SDO shows it at 30.4 nm: its chromosphere's network, active regions,
// filaments and coronal holes. Baked into a map of the surface (shaders/sdo_bake.frag,
// a few rows a frame), which the ray tracer reads: drawn inline, its registers would slow
// every ray. Needs PI, shaders/noise.glsl, and uRegionCount, uRegions[8], uGranuleFreq.

// The chromospheric network as ultraviolet and H-alpha images show it: the supergranules'
// edges, where the magnetic field gathers, a broad bright mesh broken into knots.
float filmNetwork(vec3 n, float freq, float fp) {
    float fade = 1.0 - smoothstep(0.2, 0.5, fp * freq);
    if (fade <= 0.0) return 0.0;
    vec3 p = n * freq;
    vec3 warp = vec3(gradientNoise(p * 0.7 + 5.0), gradientNoise(p * 0.7 + 9.0), gradientNoise(p * 0.7 + 2.0));
    vec3 w = worley(p + 0.9 * warp + 0.25 * vec3(gradientNoise(p * 3.1), gradientNoise(p * 3.1 + 4.0), 0.0), 0.0);
    float edge = 1.0 - smoothstep(0.0, 0.5, w.y - w.x);
    float knots = 0.3 + 0.7 * smoothstep(-0.3, 0.5, gradientNoise(p * 3.0 + 1.7));
    return fade * edge * edge * knots;
}

// Quiescent prominences (latitude, longitude in the rotating frame, size in radians, and
// 1 for a hedgerow sheet, 0 for a tall arch): two on the limbs of the first view.
const int kQuiescentProminences = 4;
const vec4 kQuiescent[kQuiescentProminences] = vec4[](vec4(0.42, -1.5, 0.32, 0.0), vec4(0.9, 1.62, 0.17, 1.0),
                                                      vec4(-0.95, 2.4, 0.15, 1.0), vec4(-0.75, -1.75, 0.13, 1.0));
// --- The Sun as SDO shows it ---------------------------------------------------------------
//
// Films and NASA's pictures show the Sun through the Solar Dynamics Observatory's 30.4 nm
// channel (He II: gas at ~80,000 K in the chromosphere and transition region), in its
// black-red-orange-yellow-white color table. There the supergranule network is bright and
// mottled by spicules and fibrils; active regions blaze with plage threaded by the bases of
// their loops; filaments (prominences seen from above) are dark threads; coronal holes over
// the poles are dim; and the limb brightens, because the line of sight runs through more
// of the thin emitting layer. Intensity 1 is the quiet network.
const int kExtraFilaments = 7;


// The local frame of a feature centered on c (rotating frame): b along it, east tilted
// toward the equator by `joy` times the latitude, and m across it.
void featureFrame(vec3 c, float joy, out vec3 b, out vec3 m) {
    vec3 east = cross(vec3(0.0, 0.0, 1.0), c);
    east = dot(east, east) > 1e-8 ? normalize(east) : vec3(1.0, 0.0, 0.0);
    vec3 north = cross(c, east);
    float tilt = joy * asin(clamp(c.z, -1.0, 1.0));
    b = east * cos(tilt) - north * sin(tilt);
    m = cross(c, b);
}

// Filaments: the quiescent prominences seen from above, and others along the filament
// channels between the active belts and the poles: wavy dark threads, tapered at the ends.
float filamentDarkness(vec3 n, float fp) {
    float dark = 0.0;
    for (int k = 0; k < kQuiescentProminences + kExtraFilaments; ++k) {
        vec3 c;
        float len;
        if (k < kQuiescentProminences) {
            vec4 q = kQuiescent[k];
            c = vec3(cos(q.x) * cos(q.y), cos(q.x) * sin(q.y), sin(q.x));
            len = 0.8 * q.z;
        } else {
            vec3 h = hash33(vec3(float(k), 17.0, 3.0));
            float lat = (h.x < 0.5 ? -1.0 : 1.0) * (0.35 + 0.6 * h.y);
            float lon = 2.0 * PI * h.z;
            c = vec3(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat));
            len = 0.1 + 0.15 * fract(h.y * 7.3);
        }
        if (dot(n, c) < cos(len + 0.05)) continue;
        vec3 b, m;
        featureFrame(c, 0.25, b, m);
        float along = atan(dot(n, b), dot(n, c));
        float u = abs(along) / len;
        if (u > 1.0) continue;
        float across = dot(n, m) + 0.012 * gradientNoise(vec3(along * 25.0, float(k), 1.0));
        float width = (0.008 + 0.006 * gradientNoise(vec3(along * 14.0, float(k), 4.0))) * (1.0 - u * u * u * u);
        width = max(width, 0.0);
        float w = max(width, fp);
        dark = max(dark, exp(-across * across / (w * w)) * (width / w));
    }
    return dark;
}

// Brightness at 30.4 nm at surface point n (rotating frame), fp the pixel's size there,
// time the chromosphere's evolution (the limb brightening is the ray tracer's).
float sdoIntensity(vec3 n, float fp, float time) {
    float freq = uGranuleFreq / 25.0;  // supergranules
    float mesh = filmNetwork(n, freq, fp);
    float mottle = fbmFiltered(n * freq * 2.5 + vec3(31.0, 0.0, time), 5, fp * freq * 2.5);
    float fibrils = ridgedFiltered(n * freq * 1.1 + vec3(7.0, time * 0.5, 0.0), 4, fp * freq * 1.1);
    // Active regions: plage in streaks along the bipole (the bright feet of its loops).
    float plage = 0.0;
    for (int k = 0; k < uRegionCount; ++k) {
        vec3 c = uRegions[k].xyz;
        float r = uRegions[k].w;
        float d = sqrt(max(2.0 - 2.0 * dot(n, c), 0.0)) / r;
        if (d > 2.5) continue;
        vec3 b, m;
        featureFrame(c, 0.5, b, m);
        vec2 local = vec2(dot(n - c, b), dot(n - c, m)) / r;
        float streaks = fbmFiltered(vec3(local.x * 3.0, local.y * 9.0, float(k) * 3.7 + time * 0.3), 4, fp / r * 9.0);
        float knots = ridgedFiltered(vec3(local * 5.0, float(k) * 1.3 + time * 0.2), 3, fp / r * 5.0);
        plage = max(plage, exp(-d * d * 1.6) * (0.7 + 0.7 * streaks + 0.6 * knots) + 0.3 * exp(-d * d * 0.5));
    }
    float hole = smoothstep(0.72, 0.92, abs(n.z)) * (0.6 + 0.4 * gradientNoise(n * 3.0 + 5.0));
    // Enhanced network: decaying old active regions, patches of brighter, denser network.
    vec3 q = n * 3.5 + vec3(0.0, 0.0, time * 0.05);
    vec3 warp = vec3(gradientNoise(q + 2.0), gradientNoise(q + 9.0), gradientNoise(q + 4.0));
    float enhanced = smoothstep(0.05, 0.55, fbmFiltered(q + 1.2 * warp, 4, fp * 3.5));
    float I = 0.55 + 0.32 * mesh * (1.0 + 1.2 * enhanced) + 0.32 * mottle + 0.42 * fibrils + 0.35 * enhanced +
              1.25 * plage;
    I *= 1.0 - 0.45 * hole;
    I *= 1.0 - 0.75 * filamentDarkness(n, fp);
    return I;
}
