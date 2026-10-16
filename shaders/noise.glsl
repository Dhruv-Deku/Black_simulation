// Hashing and value noise.
// hash13 and hash43 are "Hash without Sine" by David Hoskins (MIT license).

float hash13(vec3 p3) {
    p3 = fract(p3 * 0.1031);
    p3 += dot(p3, p3.zyx + 31.32);
    return fract((p3.x + p3.y) * p3.z);
}

vec4 hash43(vec3 p) {
    vec4 p4 = fract(vec4(p.xyzx) * vec4(0.1031, 0.1030, 0.0973, 0.1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    return fract((p4.xxyz + p4.yzzw) * p4.zywx);
}

float valueNoise(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    vec3 u = f * f * f * (f * (f * 6.0 - 15.0) + 10.0);
    float n000 = hash13(i);
    float n100 = hash13(i + vec3(1, 0, 0));
    float n010 = hash13(i + vec3(0, 1, 0));
    float n110 = hash13(i + vec3(1, 1, 0));
    float n001 = hash13(i + vec3(0, 0, 1));
    float n101 = hash13(i + vec3(1, 0, 1));
    float n011 = hash13(i + vec3(0, 1, 1));
    float n111 = hash13(i + vec3(1, 1, 1));
    return mix(mix(mix(n000, n100, u.x), mix(n010, n110, u.x), u.y),
               mix(mix(n001, n101, u.x), mix(n011, n111, u.x), u.y), u.z);
}

// Five octaves, normalized to roughly [0, 1].
float fbm(vec3 p) {
    float sum = 0.0;
    float amplitude = 0.5;
    for (int i = 0; i < 5; ++i) {
        sum += amplitude * valueNoise(p);
        p = p * 2.07 + vec3(13.7, 7.3, 3.1);
        amplitude *= 0.5;
    }
    return sum / 0.96875;
}

// --- Smoother noise for surfaces ---------------------------------------------------
// hash33 is also David Hoskins' "Hash without Sine" (MIT license).

vec3 hash33(vec3 p3) {
    p3 = fract(p3 * vec3(0.1031, 0.1030, 0.0973));
    p3 += dot(p3, p3.yxz + 33.33);
    return fract((p3.xxy + p3.yxx) * p3.zyx);
}

// Gradient (Perlin-style) noise in about [-1, 1]: random gradients at the lattice
// points, blended with a quintic fade. Unlike value noise it shows no grid.
float gradientNoise(vec3 p) {
    vec3 i = floor(p);
    vec3 f = p - i;
    vec3 u = f * f * f * (f * (f * 6.0 - 15.0) + 10.0);
    float n000 = dot(hash33(i) * 2.0 - 1.0, f);
    float n100 = dot(hash33(i + vec3(1, 0, 0)) * 2.0 - 1.0, f - vec3(1, 0, 0));
    float n010 = dot(hash33(i + vec3(0, 1, 0)) * 2.0 - 1.0, f - vec3(0, 1, 0));
    float n110 = dot(hash33(i + vec3(1, 1, 0)) * 2.0 - 1.0, f - vec3(1, 1, 0));
    float n001 = dot(hash33(i + vec3(0, 0, 1)) * 2.0 - 1.0, f - vec3(0, 0, 1));
    float n101 = dot(hash33(i + vec3(1, 0, 1)) * 2.0 - 1.0, f - vec3(1, 0, 1));
    float n011 = dot(hash33(i + vec3(0, 1, 1)) * 2.0 - 1.0, f - vec3(0, 1, 1));
    float n111 = dot(hash33(i + vec3(1, 1, 1)) * 2.0 - 1.0, f - vec3(1, 1, 1));
    return 1.6 * mix(mix(mix(n000, n100, u.x), mix(n010, n110, u.x), u.y),
                     mix(mix(n001, n101, u.x), mix(n011, n111, u.x), u.y), u.z);
}

// Filtered fractal noise in about [-1, 1]. `footprint` is the size of a pixel in the
// units of p: octaves finer than about two pixels fade out instead of aliasing into
// sparkle (their average is zero), so detail appears as the camera closes in.
float fbmFiltered(vec3 p, int octaves, float footprint) {
    float sum = 0.0;
    float amplitude = 0.5;
    float wavelength = 1.0;
    for (int i = 0; i < octaves; ++i) {
        float fade = 1.0 - smoothstep(0.25, 0.5, footprint / wavelength);
        if (fade <= 0.0) break;
        sum += amplitude * fade * gradientNoise(p);
        p = p * 2.03 + vec3(17.1, 9.7, 4.3);
        amplitude *= 0.5;
        wavelength *= 0.4926;
    }
    return sum;
}

// Ridged noise for mountain chains: sharp crests where gradient noise crosses zero.
float ridgedFiltered(vec3 p, int octaves, float footprint) {
    float sum = 0.0;
    float amplitude = 0.5;
    float wavelength = 1.0;
    float weight = 1.0;
    for (int i = 0; i < octaves; ++i) {
        float fade = 1.0 - smoothstep(0.25, 0.5, footprint / wavelength);
        if (fade <= 0.0) break;
        float n = 1.0 - abs(gradientNoise(p));
        n *= n * weight;
        weight = clamp(n * 1.5, 0.0, 1.0);
        sum += amplitude * fade * n;
        p = p * 2.07 + vec3(5.3, 11.9, 2.1);
        amplitude *= 0.5;
        wavelength *= 0.4831;
    }
    return sum;
}

// Cellular (Worley) noise: distances to the nearest and second-nearest of one random
// point per unit cell, and a random number for the nearest point's cell. The points
// drift on small loops as `phase` advances (one loop per unit), so a pattern of cells
// such as solar granulation keeps changing without ever jumping.
vec3 worley(vec3 p, float phase) {
    vec3 i = floor(p);
    vec3 f = p - i;
    float f1 = 8.0, f2 = 8.0, id = 0.0;
    for (int z = -1; z <= 1; ++z) {
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                vec3 o = vec3(x, y, z);
                vec3 h = hash33(i + o);
                vec3 point = o + 0.5 + 0.38 * sin(6.2831853 * (h + phase * (0.6 + 0.8 * h.zxy)));
                vec3 d = point - f;
                float dd = dot(d, d);
                if (dd < f1) {
                    f2 = f1;
                    f1 = dd;
                    id = h.x;
                } else if (dd < f2) {
                    f2 = dd;
                }
            }
        }
    }
    return vec3(sqrt(f1), sqrt(f2), id);
}
