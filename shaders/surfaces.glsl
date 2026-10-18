// Materials of stars and planets: what their surfaces look like, from the physics
// that shapes them (src/surface.h is the tested CPU side that sets the parameters).
// Used for the central object and for orbiting bodies alike.
//
// Points on a surface are unit vectors n in the body's rotating frame (z along the
// spin axis). `fp` is the size of a pixel on the surface in radians: detail finer than
// a pixel or two fades out (filtered noise), so the same surface stays calm seen
// whole and resolves into convection cells or coastlines as the camera closes in.

// --- Photospheres ------------------------------------------------------------------

// Eddington-Barbier for a grey atmosphere: light leaving at mu = cos(angle from the
// normal) comes from optical depth ~mu, where T^4 = (3/4) Teff^4 (tau + 2/3). So the
// limb shows cooler, higher layers: darker, and redder (the Planck curve falls faster
// in the blue). Integrated over the spectrum this is exactly Eddington's law
// I(mu) / I(1) = (2 + 3 mu) / 5; within ~10% of the measured solar limb darkening at
// every visible wavelength.
float limbTemperature(float mu) {
    return pow(0.75 * (clamp(mu, 0.0, 1.0) + 2.0 / 3.0), 0.25);
}
// At disk center, T(1) = (5/4)^(1/4) Teff.
const float kCenterTemperature = 1.05737126;

// One layer of convection cells: bright where hot gas wells up, dark in the lanes where
// it sinks back. 0..~1.2.
float convectionCells(vec3 p, float phase, float laneWidth) {
    // Bend the cell walls so the cells are rounded polygons, not crystal facets.
    vec3 warp = vec3(gradientNoise(p * 0.45 + 3.1), gradientNoise(p * 0.45 + 7.7), gradientNoise(p * 0.45 + 1.3));
    vec3 w = worley(p + (0.3 + 0.8 * laneWidth) * warp, phase);
    // Measured as ratios, which do not care how far a cell's point lies off the
    // surface (the cells are 3D, sliced by the surface): 0 on the wall between two
    // cells, rising toward each cell's middle.
    float ratio = w.x / max(w.y, 1e-4);
    float inside = 1.0 - ratio;
    float lane = smoothstep(0.0, laneWidth * 1.4, inside);
    if (laneWidth > 0.3) {
        // Giants: broad, fuzzy lanes, and many fade out, merging cells into larger,
        // irregular bright regions.
        lane = mix(1.0, sqrt(lane), smoothstep(-0.15, 0.25, gradientNoise(p * 0.6 + 31.0)));
    }
    // Each cell rises and fades over its life; the biggest split, with a dark center
    // where the gas stalls ("exploding granules").
    float life = 0.5 + 0.5 * sin(6.2831853 * (w.z * 3.7 + phase * (0.6 + 0.5 * w.z)));
    float dome = 0.75 + 0.25 * inside;
    float split = step(0.85, fract(w.z * 13.1)) * smoothstep(0.55, 0.9, inside);
    return lane * dome * (0.65 + 0.55 * life) * (1.0 - 0.45 * split);
}

// Granulation: the tops of convection cells. Returns the temperature contrast in units
// of its rms (about zero mean). freq: cells per radian (10 pressure scale heights each,
// src/surface.cpp); phase: the cells' evolution, one loop per lifetime; laneWidth: the
// lanes' width relative to the cells (narrow on dwarfs; broad on giants, whose few
// huge cells also carry smaller ones on their backs).
float granulation(vec3 n, float freq, float phase, float fp, float laneWidth) {
    // Cells under ~3 pixels across would only add noise: they fade into their average.
    float fade = 1.0 - smoothstep(0.12, 0.3, fp * freq);
    if (fade <= 0.0) return 0.0;
    vec3 p = n * freq;
    float g = convectionCells(p, phase, laneWidth);
    if (laneWidth > 0.3) {
        // Giants: the huge cells are blotchy, and carry faint small ones on their backs.
        g *= 0.65 + 0.7 * (0.5 + fbmFiltered(p * 1.1 + 4.0, 4, fp * freq * 1.1));
        float small = 1.0 - smoothstep(0.12, 0.3, fp * freq * 5.0);
        if (small > 0.0) g *= 1.0 + small * 0.12 * (convectionCells(p * 5.0 + 17.0, phase * 3.0, 0.25) - 0.6);
    }
    // Faint texture inside the cells (fine structure, when resolved).
    g += 0.1 * fbmFiltered(p * 3.0, 3, fp * freq * 3.0);
    return fade * (g - 0.62) / 0.28;
}

// The supergranulation network: cells ~25 times larger, whose edges collect magnetic
// field. Bright (faculae) toward the limb. Returns 0..1, 1 on the network.
float supergranuleNetwork(vec3 n, float freq, float fp) {
    float fade = 1.0 - smoothstep(0.15, 0.4, fp * freq);
    if (fade <= 0.0) return 0.0;
    vec3 p = n * freq;
    vec3 warp = vec3(gradientNoise(p * 0.6 + 5.0), gradientNoise(p * 0.6 + 9.0), gradientNoise(p * 0.6 + 2.0));
    vec3 w = worley(p + 0.4 * warp, 0.0);
    float edge = 1.0 - smoothstep(0.0, 0.12, w.y - w.x);
    return fade * edge * (0.6 + 0.4 * gradientNoise(p * 4.0));
}

// Starspots: where strong magnetic fields choke convection the surface is cooler. A
// dark umbra in a lighter penumbra of radial filaments; small spots (pores) have no
// penumbra. Spots come from a list (positions and sizes from src/surface.cpp, in the
// butterfly-diagram latitudes and bipolar groups of a real star: uSpots). Returns the
// temperature factor, and how much of the point is penumbra.
float spotTemperature(vec3 n, float drop, float fp, out float penumbraOut) {
    float factor = 1.0;
    penumbraOut = 0.0;
    for (int k = 0; k < uSpotCount; ++k) {
        float c = dot(n, uSpots[k].xyz);
        float radius = uSpots[k].w;
        if (c < cos(radius * 1.4)) continue;
        float d = sqrt(max(2.0 - 2.0 * c, 0.0));  // chord: the angle, for small spots
        // Ragged edges.
        float edge = 1.0 + 0.18 * gradientNoise(n * (6.0 / radius) + float(k) * 7.3);
        float x = d / (radius * edge);
        if (x >= 1.15) continue;
        float soft = max(fp / radius, 0.03);  // never sharper than a pixel
        float umbraR = uSpotUmbra[k];
        float inSpot = 1.0 - smoothstep(1.0 - soft, 1.0 + soft, x);
        float inUmbra = 1.0 - smoothstep(umbraR * (1.0 - soft), umbraR * (1.0 + soft) + 0.02, x);
        // Penumbral filaments: bright and dark streaks running out from the umbra.
        vec3 axis = uSpots[k].xyz;
        vec3 e1 = normalize(cross(axis, abs(axis.z) < 0.9 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0)));
        vec3 e2 = cross(axis, e1);
        float az = atan(dot(n, e2), dot(n, e1));
        float streaks = gradientNoise(vec3(cos(az) * 9.0, sin(az) * 9.0, x * 1.5 + float(k)));
        float filament = (1.0 - smoothstep(0.3, 0.6, fp / (radius * 0.15))) * 0.35 * streaks;
        float penumbra = 1.0 - drop * (0.17 + 0.08 * filament);
        float umbra = 1.0 - drop * (1.0 + 0.05 * gradientNoise(n * (40.0 / radius)));
        float spotFactor = mix(1.0, mix(penumbra, umbra, inUmbra), inSpot);
        factor = min(factor, spotFactor);
        penumbraOut = max(penumbraOut, inSpot * (1.0 - inUmbra));
    }
    return factor;
}

// Faculae and plage: bright magnetic patches around active regions, contrast rising
// toward the limb (they are hot walls of small flux tubes, seen from the side).
float faculae(vec3 n, float network) {
    float f = 0.0;
    for (int k = 0; k < uRegionCount; ++k) {
        float c = dot(n, uRegions[k].xyz);
        float r = uRegions[k].w;
        float d = sqrt(max(2.0 - 2.0 * c, 0.0)) / r;
        f = max(f, exp(-d * d * 1.5));
    }
    return f * (0.45 + 0.55 * network);
}

// A brown dwarf's weather: too cool to shine like a star, its atmosphere forms clouds
// of silicate and iron in latitude bands, like Jupiter's, which we see in thermal
// glow: thin clouds let out light from deeper, hotter layers. Returns a temperature
// factor (observed brightness swings of ~10% at 1,300 K mean ~2-3% in temperature).
float brownDwarfWeather(vec3 n, float fp, float seed) {
    float lat = asin(clamp(n.z, -1.0, 1.0));
    float turb = fbmFiltered(n * vec3(3.0, 3.0, 9.0) + seed, 6, fp * 9.0);
    float bands = sin(lat * 14.0 + 2.2 * turb);
    float patches = fbmFiltered(n * vec3(6.0, 6.0, 14.0) + seed * 1.7 + 13.0, 5, fp * 14.0);
    // A large storm or two, like Jupiter's spots.
    vec3 stormDir = normalize(vec3(cos(seed * 3.1), sin(seed * 3.1), 0.25));
    float storm = exp(-pow(length((n - stormDir) * vec3(1.0, 1.0, 2.3)) / 0.16, 2.0));
    return 1.0 + 0.022 * bands + 0.018 * patches - 0.03 * storm;
}

// --- Earth-like planets ------------------------------------------------------------

// A rotation of p about the unit axis k by angle a (Rodrigues).
vec3 rotateAbout(vec3 p, vec3 k, float a) {
    float c = cos(a), s = sin(a);
    return p * c + cross(k, p) * s + k * dot(k, p) * (1.0 - c);
}

// The nearest feature point of cellular noise (not just its distance), for placing
// storms.
vec3 nearestCellPoint(vec3 p) {
    vec3 i = floor(p);
    vec3 best = vec3(0.0);
    float bestD = 1e9;
    for (int z = -1; z <= 1; ++z) {
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                vec3 c = i + vec3(x, y, z);
                vec3 point = c + 0.15 + 0.7 * hash33(c);
                float d = dot(point - p, point - p);
                if (d < bestD) {
                    bestD = d;
                    best = point;
                }
            }
        }
    }
    return best;
}

// Cloud cover of an Earth-like planet, 0..1: by latitude (the rainy intertropical band
// just north of the equator, clear subtropics where the air sinks, stormy
// mid-latitudes), shaped by swirling cyclones in the mid-latitudes (counter-clockwise
// in the north, clockwise in the south).
float earthClouds(vec3 n, float fp, float seed) {
    float lat = asin(clamp(n.z, -1.0, 1.0));
    float absLat = abs(lat);
    float cq1 = (lat - 0.09) / 0.07, cq2 = (absLat - 0.44) / 0.11, cq3 = (absLat - 1.0) / 0.25;
    float cover = 0.48 + 0.24 * exp(-cq1 * cq1) - 0.2 * exp(-cq2 * cq2) + 0.2 * exp(-cq3 * cq3);
    // Mid-latitude cyclones wind the cloud field into commas. Storm centers: points of
    // a cellular pattern, projected onto the sphere.
    vec3 p = n;
    vec3 center = normalize(nearestCellPoint(n * 2.6 + seed * 0.3) - vec3(seed * 0.3));
    float cLat = asin(clamp(center.z, -1.0, 1.0));
    float dc = acos(clamp(dot(n, center), -1.0, 1.0));
    if (abs(cLat) > 0.6 && abs(cLat) < 1.15 && dc < 0.32) {
        // Sampling the noise at points turned by +a(d) shows it turned by -a(d): arms that
        // wind counter-clockwise inward (air spiralling into a northern low) need a < 0.
        p = rotateAbout(n, center, -2.6 * pow(1.0 - dc / 0.32, 1.5) * sign(cLat));
    }
    // The winds stretch clouds along the latitude circles, most in the jet streams.
    float stretch = 1.0 + 1.6 * smoothstep(0.4, 0.8, absLat);
    vec3 q = p * vec3(4.0, 4.0, 4.0 * stretch);
    float fq = fp * 4.0 * stretch;
    vec3 warp = vec3(fbmFiltered(q * 0.5 + 41.0, 4, fq * 0.5), fbmFiltered(q * 0.5 + 12.0, 4, fq * 0.5),
                     fbmFiltered(q * 0.5 + 77.0, 4, fq * 0.5));
    float c = fbmFiltered(q + 2.5 * warp + seed * 1.9, 12, fq);
    c += 0.25 * fbmFiltered(p * 18.0 + 3.0 * warp + 5.0, 10, fp * 18.0);
    // Fields of small cumulus cells, where resolved.
    float cellFade = 1.0 - smoothstep(0.12, 0.3, fp * 70.0);
    if (cellFade > 0.0) c += cellFade * 0.1 * (0.45 - worley(p * 70.0 + seed, 0.0).x);
    // c is roughly normal with sd ~0.2: this threshold covers about `cover` of the sky.
    // Past it the optical depth tau grows; a cloud reflects ~tau / (tau + 6) (two-stream,
    // forward-scattering droplets), so thin edges stay translucent.
    float tau = max(c - 0.42 * (0.5 - cover), 0.0) * 60.0;
    return tau / (tau + 6.0);
}

struct Ground {
    vec3 albedo;    // diffuse reflectance (linear sRGB), hill shading included
    float water;    // open water: mirror-like sun glint
    float clouds;   // cloud cover, 0..1 (drawn over the ground)
    float lights;   // city lights at night
};

// The land and sea of an Earth-like planet. Continents from warped fractal noise, set
// so that ~30% is land; mountain chains from ridged noise; climate from latitude: wet
// tropics, desert belts near 25 degrees (where the Hadley cells' air sinks), forests in
// the wet mid-latitudes, tundra and ice toward the poles and up high. Reflectances are
// those of the real surfaces (dark forest ~0.05-0.1, desert ~0.35, snow ~0.8, deep
// ocean ~0.03). sun: direction to the sun in the same frame as n (zero: no hill
// shading).
Ground earthGround(vec3 n, float fp, float seed, vec3 sun) {
    Ground gr;
    float lat = asin(clamp(n.z, -1.0, 1.0));
    float absLat = abs(lat);
    vec3 q = n * 1.3 + seed;
    vec3 warp = vec3(gradientNoise(q + 1.7), gradientNoise(q + 9.2), gradientNoise(q + 4.4));
    // Octaves down to ~100 m: they only cost anything when the camera is close.
    float height = fbmFiltered(n * 1.7 + 0.6 * warp + seed, 16, fp * 1.7);
    float relief = ridgedFiltered(n * 5.0 + 0.7 * warp + seed * 2.0, 14, fp * 5.0);
    const float sea = 0.06;
    float elevation = height - sea + 0.12 * relief * smoothstep(-0.02, 0.1, height - sea);
    float coastBlur = max(fp * 0.6, 0.002);
    float land = smoothstep(-coastBlur, coastBlur, elevation);

    // Climate.
    float sq1 = absLat / 0.17, sq2 = (absLat - 0.44) / 0.11, sq3 = (absLat - 0.95) / 0.22;
    float climate = 0.35 * exp(-sq1 * sq1) - 0.45 * exp(-sq2 * sq2) + 0.15 * exp(-sq3 * sq3);
    float wet = clamp(0.55 + climate + 0.5 * fbmFiltered(n * 3.5 + 31.0, 6, fp * 3.5) -
                          0.15 * smoothstep(0.1, 0.4, elevation),  // dry continental interiors
                      0.0, 1.0);
    // Colder toward the poles, and high up (only real mountain ranges matter).
    float temperature = 1.0 - absLat / 1.25 - 0.55 * smoothstep(0.15, 0.5, elevation) +
                        0.06 * gradientNoise(n * 6.0 + 3.0);

    float tint = fbmFiltered(n * 11.0 + 51.0, 13, fp * 11.0);
    vec3 desert = mix(vec3(0.46, 0.31, 0.16), vec3(0.62, 0.48, 0.31), 0.5 + 0.5 * tint);
    vec3 steppe = vec3(0.23, 0.20, 0.11);
    vec3 grass = vec3(0.09, 0.12, 0.045);
    vec3 forest = vec3(0.032, 0.062, 0.022);
    vec3 tundra = vec3(0.12, 0.11, 0.085);
    vec3 rock = vec3(0.17, 0.15, 0.13);
    vec3 snow = vec3(0.78, 0.81, 0.86);
    vec3 ground = mix(desert, steppe, smoothstep(0.18, 0.36, wet));
    ground = mix(ground, grass, smoothstep(0.36, 0.52, wet));
    ground = mix(ground, forest, smoothstep(0.52, 0.72, wet));
    ground = mix(ground, tundra, (1.0 - smoothstep(0.16, 0.32, temperature)));
    ground = mix(ground, rock, smoothstep(0.25, 0.45, relief) * smoothstep(0.2, 0.45, elevation));
    ground *= 0.85 + 0.3 * (0.5 + 0.5 * tint);
    // Fields, forest patches, river valleys: texture that only shows close up.
    ground *= 1.0 + 0.9 * fbmFiltered(n * 900.0 + 7.0, 9, fp * 900.0);
    float ice = (1.0 - smoothstep(0.02, 0.1, temperature + 0.05 * tint));
    ground = mix(ground, snow, ice);
    // Hill shading: slopes facing the sun are brighter, those facing away darker (from
    // the mountains' noise one step toward the sun; exaggerated, since the real ~1 km
    // relief over 6,371 km would hardly show).
    float hills = 1.0;
    if (dot(sun, sun) > 0.0 && land > 0.0) {
        vec3 t = sun - n * dot(sun, n);
        float tl = length(t);
        if (tl > 1e-4) {
            float e = max(fp, 3e-5);
            float ahead = ridgedFiltered((n + t / tl * e) * 5.0 + 0.7 * warp + seed * 2.0, 14, fp * 5.0);
            float slope = (relief - ahead) / e * smoothstep(0.0, 0.08, elevation);
            hills = clamp(1.0 + 0.004 * slope, 0.5, 1.6);
        }
    }
    ground *= hills;

    // Ocean: dark blue, turquoise over shallow shelves, sea ice near the poles.
    float shelf = smoothstep(-0.06, 0.0, elevation);
    vec3 ocean = mix(vec3(0.004, 0.016, 0.042), vec3(0.02, 0.075, 0.085), shelf * shelf);
    float seaIce = smoothstep(1.17, 1.25, absLat + 0.04 * gradientNoise(n * 9.0 + 7.0));
    ocean = mix(ocean, vec3(0.65, 0.7, 0.76), seaIce);

    gr.albedo = mix(ocean, ground, land);
    gr.water = (1.0 - land) * (1.0 - seaIce);

    // Night lights: cities with bright cores and scattered towns, in mild, wet-enough
    // lowlands. Too small to see, they blend into their average glow.
    float habitable = land * smoothstep(0.25, 0.45, wet) * smoothstep(0.25, 0.45, temperature) *
                      (1.0 - smoothstep(0.15, 0.3, elevation)) * (1.0 - smoothstep(0.0, 0.75, ice));
    float lights = 0.0;
    if (habitable > 0.0) {
        float density = smoothstep(-0.05, 0.3, fbmFiltered(n * 7.0 + 71.0, 5, fp * 7.0));
        vec3 c1 = worley(n * 45.0 + 3.0, 0.0);
        float city = step(c1.z, 0.45) * exp(-c1.x * c1.x / 0.02);
        city = mix(0.028, city, 1.0 - smoothstep(0.1, 0.3, fp * 45.0));
        vec3 c2 = worley(n * 180.0 + 9.0, 0.0);
        float town = step(c2.z, 0.3) * exp(-c2.x * c2.x / 0.03);
        town = mix(0.028, town, 1.0 - smoothstep(0.1, 0.3, fp * 180.0));
        lights = density * (3.0 * city + town);
    }
    gr.lights = habitable * lights;
    gr.clouds = earthClouds(n, fp, seed);
    return gr;
}
// --- Gas giants ----------------------------------------------------------------------

// Jupiter's face. Its belts (dark, sinking air, reddish-brown) and zones (bright,
// rising, ammonia-ice white) at their real latitudes; the shear of the jets between
// them rolls the boundaries into eddies; the Great Red Spot at 22 degrees south (an
// anticyclone ~16,000 by 11,000 km, turning counter-clockwise), white ovals, the
// blue-grey "hot spots" and festoons at the edge of the North Equatorial Belt, and
// the mottled grey-blue polar regions seen by Juno.
vec3 jupiterAlbedo(vec3 n, float fp, float seed) {
    float lat = asin(clamp(n.z, -1.0, 1.0));
    float lon = atan(n.y, n.x);
    // The Great Red Spot: an ellipse 12 by 8 degrees; inside, the flow spirals.
    const float grsLat = -0.384;
    float grsLon = seed * 2.0 - 2.4;
    vec2 grs = vec2(lon - grsLon, lat - grsLat);
    grs.x = mod(grs.x + 3.14159265, 6.2831853) - 3.14159265;
    vec2 e = grs / vec2(0.105, 0.072);
    float grsR = length(e);
    vec3 pn = n;
    if (grsR < 1.8) {
        // Counter-clockwise seen from outside: an anticyclone in the south.
        vec3 axis = vec3(cos(grsLat) * cos(grsLon), cos(grsLat) * sin(grsLon), sin(grsLat));
        pn = rotateAbout(n, axis, 2.4 * exp(-grsR * grsR * 0.9));
    }
    // Shear: eddies roll up between the jets. The noise is stretched along longitude and
    // warped twice through itself (Quilez's domain warping), which folds it into the
    // swirls and filaments of a sheared fluid.
    vec3 s = pn * vec3(2.5, 2.5, 9.0) + seed;
    float fs = fp * 9.0;
    vec3 q = vec3(fbmFiltered(s, 5, fs), fbmFiltered(s + vec3(5.2, 1.3, 2.8), 5, fs),
                  fbmFiltered(s + vec3(1.7, 9.2, 4.1), 5, fs));
    vec3 w2 = vec3(fbmFiltered(s + 4.0 * q + vec3(1.7, 9.2, 0.3), 5, fs),
                   fbmFiltered(s + 4.0 * q + vec3(8.3, 2.8, 5.1), 5, fs), 0.0);
    float turb = fbmFiltered(s + 4.0 * w2, 7, fs);  // about +-0.4
    float fine = fbmFiltered(pn * vec3(12.0, 12.0, 40.0) + 3.0 * w2 + 9.0, 7, fp * 40.0);
    float l = lat + 0.03 * turb + 0.008 * fine;  // wavy band edges
    float deg = l * 57.29578;

    // Albedo by latitude (planetographic degrees): zones bright, belts brown.
    vec3 zone = vec3(0.66, 0.61, 0.50);
    vec3 belt = vec3(0.44, 0.24, 0.12);
    vec3 col = zone;
    // Each belt: center, half-width, darkness.
    col = mix(col, belt, 0.9 * (1.0 - smoothstep(3.0, 5.5, abs(deg - 12.5))));    // North Equatorial Belt
    col = mix(col, belt, 0.85 * (1.0 - smoothstep(4.0, 6.5, abs(deg + 13.5))));   // South Equatorial Belt
    col = mix(col, belt * 1.15, 0.7 * (1.0 - smoothstep(1.8, 3.5, abs(deg - 27.5))));  // North Temperate Belt
    col = mix(col, belt * 1.15, 0.6 * (1.0 - smoothstep(1.4, 3.0, abs(deg + 30.0))));  // South Temperate Belt
    col = mix(col, belt * 1.3, 0.45 * (1.0 - smoothstep(1.0, 2.5, abs(deg - 37.0))));   // North North Temperate Belt
    col = mix(col, belt * 1.3, 0.4 * (1.0 - smoothstep(1.0, 2.5, abs(deg + 41.0))));
    col = mix(col, vec3(0.80, 0.74, 0.60), 0.4 * (1.0 - smoothstep(0.0, 4.0, abs(deg))));  // Equatorial Zone
    // Orange haze in the Equatorial Zone some years.
    col *= mix(vec3(1.0), vec3(1.0, 0.9, 0.76), 0.5 * (1.0 - smoothstep(2.0, 6.0, abs(deg))));
    // Polar regions: grey-blue, mottled with cyclones.
    float polar = smoothstep(45.0, 62.0, abs(deg));
    float mottle = fbmFiltered(n * 10.0 + 3.0 * q + 21.0, 7, fp * 10.0);
    col = mix(col, vec3(0.30, 0.30, 0.32) * (0.8 + 0.9 * mottle), polar);
    // Swirls and streaks everywhere: the belts' dark filaments, the zones' bright puffs.
    col *= 1.0 + 0.6 * turb + 0.15 * fine;
    // Festoons and blue-grey hot spots on the south edge of the North Equatorial Belt.
    float festoonBand = (1.0 - smoothstep(0.5, 3.0, abs(deg - 7.0)));
    float hotspots = smoothstep(0.35, 0.7, gradientNoise(vec3(lon * 4.5, deg * 0.3, seed)));
    col = mix(col, vec3(0.30, 0.34, 0.40), festoonBand * hotspots * 0.75);

    // The Great Red Spot: brick-red core, paler collar, the "hollow" it carves in the
    // South Equatorial Belt.
    float core = 1.0 - smoothstep(0.55, 1.0, grsR);
    float collar = smoothstep(0.75, 1.0, grsR) * (1.0 - smoothstep(1.0, 1.35, grsR));
    col = mix(col, vec3(0.74, 0.68, 0.58), collar * 0.8);
    col = mix(col, vec3(0.66, 0.24, 0.10) * (0.85 + 0.5 * turb), core * 0.95);
    // White ovals in the South Temperate Belt.
    for (int k = 0; k < 3; ++k) {
        vec2 o = vec2(lon - (seed + 1.4 + 1.9 * float(k)), lat + 0.58);
        o.x = mod(o.x + 3.14159265, 6.2831853) - 3.14159265;
        float d = length(o / vec2(0.04, 0.025));
        col = mix(col, vec3(0.88, 0.86, 0.80), (1.0 - smoothstep(0.7, 1.0, d)) * 0.9);
    }
    return col;
}

// Minnaert's law for a hazy, cloudy planet: reflected radiance ~ mu0^k mu^(k-1)
// (k = 1 is Lambert). Jupiter's clouds in visible light have k ~ 0.9: the limb darkens
// a little less steeply than for a Lambert sphere, seen at full phase.
float minnaert(float mu0, float mu, float k) {
    return pow(max(mu0, 0.0), k) * pow(max(mu, 0.05), k - 1.0);
}

// Sun glint off water: the GGX microfacet model (wind-roughened sea, slope variance
// ~0.03 for a 5 m/s wind, Cox and Munk) with Fresnel reflection from water
// (n = 1.33, 2% at normal incidence). Returns the reflected radiance per unit
// incident radiance (pi f mu0, matching albedo * mu0 for a diffuse surface).
float sunGlint(vec3 normal, vec3 toLight, vec3 toViewer, float roughness) {
    float nl = dot(normal, toLight), nv = dot(normal, toViewer);
    if (nl <= 0.0 || nv <= 0.0) return 0.0;
    vec3 h = normalize(toLight + toViewer);
    float nh = max(dot(normal, h), 0.0), vh = max(dot(toViewer, h), 0.0);
    float a2 = roughness * roughness;
    float dd = nh * nh * (a2 - 1.0) + 1.0;
    float D = a2 / (3.14159265 * dd * dd);
    float k = roughness * 0.5;
    float G = (nl / (nl * (1.0 - k) + k)) * (nv / (nv * (1.0 - k) + k));
    float F = 0.02 + 0.98 * pow(1.0 - vh, 5.0);
    return 3.14159265 * D * G * F / (4.0 * nv);
}

// --- Airless rocky and icy bodies ---------------------------------------------------
//
// Moons, Mercury, asteroids, comet nuclei, satellites, and planets nobody has seen the
// surface of (src/bodies.h, Palette). Normal albedos from spacecraft photometry: the
// Moon's highlands ~0.15 and maria ~0.07; Io's sulfur plains ~0.6, with red and black
// volcanic centers; Europa's ice ~0.65, crossed by reddish-brown lineae; Ganymede ~0.43,
// dark ancient terrain cut by bright grooved lanes; Callisto ~0.2, peppered with bright
// crater ejecta; Mercury ~0.14; Venus's cloud tops ~0.75; comet nuclei and carbonaceous
// asteroids ~0.04.
const int kPaletteMoon = 0, kPaletteIo = 1, kPaletteEuropa = 2, kPaletteGanymede = 3, kPaletteCallisto = 4,
          kPaletteMercury = 5, kPaletteVenus = 6, kPaletteComet = 7, kPaletteSatellite = 8, kPaletteRocky = 9,
          kPaletteScorched = 10;

// Craters at three sizes: dark floors (x) and bright raised rims (y), each set faded out
// once its cells shrink below a few pixels.
vec2 craters(vec3 n, float fp, float seed) {
    vec2 sum = vec2(0.0);
    float scale = 4.0;
    for (int k = 0; k < 3; ++k) {
        float fade = 1.0 - smoothstep(0.15, 0.35, fp * scale);
        if (fade <= 0.0) break;
        // Warped, so that craters are not perfect circles on a grid.
        vec3 p = n * scale + seed * (1.0 + float(k));
        vec3 w = worley(p + 0.25 * vec3(gradientNoise(p * 1.7), gradientNoise(p * 1.7 + 5.1), 0.0), 0.0);
        float present = step(0.62, w.z);  // not every cell holds a crater
        float x = w.x / (0.12 + 0.3 * fract(w.z * 7.3));
        sum += present * fade * vec2(1.0 - smoothstep(0.5, 1.0, x), exp(-(x - 1.0) * (x - 1.0) / 0.02));
        scale *= 3.1;
    }
    return sum;
}

vec3 rockyAlbedo(vec3 n, float fp, float seed, int palette) {
    if (palette == kPaletteSatellite) return vec3(0.6);  // solar panels and foil, unresolved
    float big = fbmFiltered(n * 2.0 + seed, 5, fp * 2.0);
    float mid = fbmFiltered(n * 9.0 + seed * 1.7, 5, fp * 9.0);
    if (palette == kPaletteVenus) {
        // Cloud tops: nearly featureless in visible light, with faint banding.
        float lat = asin(clamp(n.z, -1.0, 1.0));
        return vec3(0.86, 0.80, 0.62) * (1.0 + 0.06 * fbmFiltered(vec3(n.xy * 2.0, lat * 6.0) + seed, 4, fp * 6.0));
    }
    if (palette == kPaletteIo) {
        vec3 col = mix(vec3(0.80, 0.72, 0.38), vec3(0.62, 0.38, 0.17), smoothstep(-0.1, 0.35, big));
        col = mix(col, vec3(0.86, 0.85, 0.78), smoothstep(0.15, 0.4, mid) * 0.5);  // sulfur dioxide frost
        vec3 w = worley(n * 6.0 + seed, 0.0);  // volcanic centers: dark calderas in red rings
        float vent = step(0.72, w.z);
        col = mix(col, vec3(0.55, 0.22, 0.10), vent * (1.0 - smoothstep(0.2, 0.45, w.x)) * 0.8);
        return mix(col, vec3(0.06, 0.05, 0.05), vent * (1.0 - smoothstep(0.05, 0.12, w.x)));
    }
    vec2 cr = craters(n, fp, seed);
    if (palette == kPaletteEuropa) {
        vec3 col = vec3(0.76, 0.72, 0.64) * (1.0 + 0.12 * big);
        float lineae = ridgedFiltered(n * 5.0 + seed, 5, fp * 5.0);  // crests where the ice cracked
        col = mix(col, vec3(0.50, 0.33, 0.22), smoothstep(0.55, 0.8, lineae) * 0.75);
        return mix(col, vec3(0.55, 0.42, 0.33), smoothstep(0.1, 0.4, -big) * 0.5);  // chaos terrain
    }
    if (palette == kPaletteGanymede) {
        vec3 col = mix(vec3(0.27, 0.25, 0.22), vec3(0.52, 0.51, 0.49), smoothstep(-0.1, 0.15, big + 0.3 * mid));
        return col * (1.0 + 0.3 * cr.y - 0.1 * cr.x);
    }
    if (palette == kPaletteCallisto) {
        vec3 col = vec3(0.20, 0.19, 0.17) * (1.0 + 0.2 * mid);
        return mix(col, vec3(0.6, 0.6, 0.58), clamp(cr.y * 0.8, 0.0, 1.0));  // bright fresh ejecta
    }
    // Light and dark terrain: the Moon's highlands and maria (a sixth of its surface).
    vec3 light = vec3(0.17, 0.16, 0.15), dark = vec3(0.075, 0.072, 0.07);
    if (palette == kPaletteComet) {
        light = vec3(0.05);
        dark = vec3(0.035);
    } else if (palette == kPaletteMercury) {
        light = vec3(0.15, 0.14, 0.13);
        dark = vec3(0.10, 0.095, 0.09);
    } else if (palette == kPaletteRocky) {
        light = vec3(0.22, 0.16, 0.12);
        dark = vec3(0.10, 0.08, 0.07);
    } else if (palette == kPaletteScorched) {
        light = vec3(0.12, 0.10, 0.09);
        dark = vec3(0.05, 0.045, 0.045);
    }
    // At full phase craters barely show in albedo: slightly darker floors, slightly
    // brighter fresh rims.
    vec3 col = mix(dark, light, smoothstep(-0.25, -0.08, big + 0.15 * mid));
    col *= (1.0 + 0.12 * mid) * (1.0 - 0.08 * cr.x);
    return col + light * 0.2 * cr.y;
}
