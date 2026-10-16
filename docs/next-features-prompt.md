# Build prompt: turn the Kerr black hole simulator into a compact-objects simulator

## 0. What I want, in plain words

I have a working real-time simulator of a rotating (Kerr) black hole. I want the **object at
the center to be selectable**: a black hole, a neutron star (pulsar, magnetar), a white
dwarf, a normal star like the Sun, a red dwarf, a brown dwarf, a red giant, or a planet.
Gravity, light bending, surfaces, disks, tides, accretion and everything else must change
correctly with the object. I also want a **compactness slider**: keep the mass, shrink the
radius, and watch a Sun turn into a white dwarf, then a neutron star, then collapse into a
black hole, with relativistic effects switching on as it gets smaller.

On top of that, add ten more features (section 3): the Penrose process, gravitational-wave
inspirals with a chirp you can hear, an Event Horizon Telescope view, a disk that forms from
tidal debris, naked singularities (spin above 1), timeline scrubbing, sound, cinematic camera
keyframes, saving and loading scenes, and speed-ups.

All the physics and maths you need are below. Where I give a number, use it as a test.
If anything here conflicts with the code or with physics, trust the physics, fix it, and
tell me what you changed.

---

## 1. The existing project (read this before touching anything)

**Location:** `C:\Users\DHRUV\Desktop\black hole` (Windows 11, Visual Studio 2026, CMake ≥ 3.24,
OpenGL 4.6, NVIDIA RTX 5050 Laptop GPU). GitHub: `Dhruv-Deku/Black_simulation` (may need
re-creating and pushing; do not push without asking me).

**Build and test:** `.\build.ps1` (configures, builds Release, runs the tests),
`.\build.ps1 -Run` (also starts it). Executable: `build\Release\blackhole.exe`.
The exe is locked while it runs; close it before rebuilding.

**Read first:** `README.md` (features, maths, approximations), then `src/kerr.h`,
`src/trajectory.h`, `src/bodies.h`, `src/gravity.h`, `src/tidal.h`, `src/particles.h`,
`shaders/trace.frag`, `tests/test_kerr.cpp`.

### 1.1 Architecture map

| Piece | Where | What it does |
|---|---|---|
| Units | everywhere | G = c = M = 1. M is the central mass. Lengths in GM/c², times in GM/c³. |
| Coordinates | `src/kerr.*` | Ingoing Kerr–Schild (t, r, θ, φ): r, θ are Boyer–Lindquist; t = t_BL + ∫(r²+a²)/Δ dr − r; φ = φ_BL + ∫a/Δ dr. Regular across the horizon. |
| Geodesics | `src/kerr.*` | `State{t, r, θ, φ, vr, vθ}` (vr, vθ are Mino-time velocities), `Constants{E, L, Q, μ}`. `kerr::derivatives` gives d/dτ (Mino time, dλ = Σ dτ). Second-order radial/polar equations + RK4 + constraint re-projection (`enforceConstraints`) + spin-axis snap. |
| Horizon-regular factor | `kerr::inflowFactor` | F = (P + vr)/Δ, computed as K/(P − vr) on the constraint surface when P·vr < 0. |
| Cartesian chart | `kerr::toCartesian` | Kerr–Schild Cartesian: x + iy = (r + ia) e^{iφ} sinθ, z = r cosθ. |
| Observer frames | `kerr::observerFrame`, `photonFromFrame` | Gram–Schmidt frame from any 4-velocity (ZAMO frame for a hovering camera). Photons are built with unit energy in the camera frame; E may be ≤ 0 in the ergosphere. |
| Ray tracer | `shaders/trace.frag` | Backward RK4 in Mino time; disk crossing; bodies from a history texture (position, 4-velocity, shape); debris gas from 3D textures; jets; path lines (segment grid + occupancy texture); per-ray, per-body culling margins; horizon rules (a ray from inside the horizon that reaches it with P(r+) = 2r+E − aL < 0 came from the far side → black). |
| CPU ray inspector | `src/trajectory.*` | `traceRay` (same physics, doubles), `classify` (exact fate from the zeros of R(r)), `predict`, `effectivePotential`. |
| Bodies | `src/bodies.*` | Doubles. Kick–drift–kick mutual gravity (`src/gravity.*`, rest-frame simultaneity), affine tidal model (`src/tidal.*`), light-delay history (4 rows per slot in a 4096-column texture, 0.5 M per column), departed bodies still drawn until their last light arrives. |
| Gas | `src/particles.*`, `shaders/particles_update.comp`, `particles_deposit.comp`, `grid_resolve.comp` | Each parcel follows its own geodesic; particle-mesh hydro (pressure, collisional shocks, entropy T/ρ^{γ−1}, cooling ∝ T³, floor temperature, bodies' pull); full-precision physics grid + half-precision drawing grid; light-delayed "seen" copy of every parcel. |
| Camera | `src/camera.*` | Hover / Orbit / FreeFall / RideAlong; falls through the horizon; stops just above the inner horizon (`kerr::interiorLimit`). |
| UI, app | `src/ui.*`, `src/app.*` | ImGui panel, launch planner, orbit map, ray inspector (toggle with I), CLI flags, video recording, slow motion near the hole. |
| Tests | `tests/test_kerr.cpp` | 134 checks against closed-form results. |

### 1.2 Rules of the house

- Every new piece of physics gets at least one test against a closed-form or published value.
  Section 4 lists the ones I want.
- Verify visuals headlessly: `blackhole.exe --size 1280x720 --preset … --advance T --frames N --screenshot out.png`
  (add `--ui` to include the panel), and `--record out.mp4 --duration S`. Look at the images.
- Measure performance: the screenshot run prints GPU time. The laptop GPU varies ±30 % with heat;
  compare runs side by side. Budget: about 20 ms per frame at 1080p with render scale 0.75.
- Keep comments explaining the physics, no unexplained magic numbers, and match the existing
  code style. Update `README.md` (features, maths, test count, **approximations**) as you go.
- Lessons already paid for:
  - In Git Bash, heredocs break on apostrophes; write edit scripts to files.
  - `sample` is a reserved word in GLSL.
  - Gas temperature must be carried as entropy, or density noise random-walks it upward.
  - Any force between fast-moving bodies must be computed in the body's rest frame with
    rest-frame simultaneity (a moon spun out of its planet's grip otherwise).
  - Averaged gas 4-velocities must be renormalized with the metric before computing a
    Doppler factor, and the factor capped (we use g ≤ 4).
  - Rays starting exactly in the disk plane must not count as crossing it (use a strict test).

---

## 2. Main feature: the central object

### 2.1 Units and physical scale

```
length unit  L_M = G M / c²  = 1.476625 km  × (M / M_sun)
time unit    T_M = G M / c³  = 4.925491 µs  × (M / M_sun)
speed        dimensionless (× c);  energy per unit mass in units of c²
acceleration a_SI = a_geom × c² / L_M        (a_geom in units of 1/M)
compactness  C = M / R                        (R in units of M; the Sun: 2.1e-6, a neutron star: ~0.17)
```

Add a **physical units** mode: the user gives the central mass in solar masses; every readout
shows km, AU, seconds or hours as well as M.

### 2.2 Catalogue of central objects (typical values; make each a one-click preset)

| Object | M (M_sun) | R | R/M (geometric) | spin χ = cJ/(GM²) | T_surface |
|---|---|---|---|---|---|
| Earth | 3.0e-6 | 6,371 km | 1.44e9 | ~900 | (reflects light) |
| Jupiter | 9.55e-4 | 69,911 km | 4.96e7 | ~10³ | ~165 K (reflects) |
| Brown dwarf | 0.05 | ~70,000 km | ~9.5e5 | small | 1,000–2,500 K |
| Red dwarf (Proxima) | 0.12 | 107,000 km | 6.0e5 | small | 3,000 K |
| Sun | 1 | 695,700 km | 4.71e5 | 0.22 | 5,772 K |
| Red giant | 1 | 100 R_sun | 4.7e7 | tiny | ~4,000 K |
| White dwarf (typical) | 0.6 | 8,700 km | 9,800 | 0.004 (spin period 1 day) to ~10 (33 s, AE Aqr) | 10,000–25,000 K |
| Sirius B | 1.02 | 5,840 km | 3,877 | small | 25,000 K |
| Neutron star | 1.4 | 12 km | 5.8 | up to ~0.4 (716 Hz pulsar: 0.37) | ~1e6 K (hot spots 1e7 K) |
| Stellar black hole | 10 | r+ = 29.5 km (a = 0) | 2 (a = 0) | 0–0.998 | — |
| Sgr A* | 4.3e6 | r+ | | ~0.9 | L_M = 6.35e6 km, T_M = 21.2 s |
| M87* | 6.5e9 | r+ | | ~0.9 | L_M = 9.6e9 km (64 AU), T_M = 8.9 h |

Spin check for the Sun: J_sun ≈ 1.92e41 kg m²/s, GM_sun²/c = 8.81e41 → χ = 0.22.

### 2.3 The spacetime outside each object

- **Birkhoff's theorem:** outside any spherical non-rotating body the metric is exactly
  Schwarzschild. So **for the same mass, orbits and light bending outside the surface are
  identical for a black hole, a neutron star and the Sun.** Only the surface radius changes
  how close you can get. Make this the headline teaching point in the UI.
- **Rotating compact objects:** use the existing Kerr metric with a = χ outside the surface.
  This gets the mass and angular momentum exactly. The quadrupole is wrong for a neutron
  star (Q = −q χ² M³ with q ≈ 2–10 for real neutron stars versus q = 1 for Kerr), which is a
  small effect for χ ≲ 0.4. If you want better, implement the Hartle–Thorne metric as an
  option. It is not separable (no Carter constant), so it needs the general Hamiltonian
  integrator below.
- **Non-compact bodies with χ > 1 (planets):** full Kerr with a > 1 is a naked singularity and
  must not be used as their exterior. Frame dragging near them is negligible anyway
  (ω ≈ 2J/r³), so clamp χ_eff = min(χ, 0.99) for bodies with R/M > 50 and document it.
- **Mass-shedding limit:** a star cannot spin faster than about Ω_max ≈ 0.57 √(M/R³), with M and
  R of the non-rotating star (Lattimer & Prakash 2004: f_max ≈ 1045 Hz (M/M_sun)^(1/2) (10 km/R)^(3/2)).
  Limit the spin slider with it. χ = I Ω / M², with I ≈ 0.35 M R² for a neutron star
  (check: 1.4 M_sun, 12 km, 716 Hz gives χ = 0.366).
- **General Hamiltonian geodesic integrator** (only if you add a non-separable metric):

```
H = ½ g^{μν}(x) p_μ p_ν
dx^μ/dλ = ∂H/∂p_μ = g^{μν} p_ν
dp_μ/dλ = −∂H/∂x^μ = −½ (∂_μ g^{αβ}) p_α p_β
Keep H = −½ μ² (massive) or 0 (light) as a check; use adaptive RK45 (Dormand–Prince).
```

### 2.4 Surface physics (Schwarzschild form; for Kerr use `kerr::metricAt` values)

```
Static observer at radius r (a = 0):       u = (1 − 2M/r)^(−1/2) ∂_t
Gravitational redshift from the surface:   1 + z = (1 − 2M/R)^(−1/2)
Surface gravity (proper acceleration):     g = M / ( R² √(1 − 2M/R) )
Local escape speed:                        v_esc = √(2M/R)
Surface clock rate vs. infinity:           √(1 − 2M/R)
Buchdahl limit (no static star smaller):   R > 9M/4
Photon sphere visible outside the star:    only if R < 3M   (a = 0)
ISCO outside the star:                     if R < 6M         (a = 0; typical neutron star R ≈ 5.8M → yes)
```

Worked numbers (use as tests):
- Neutron star 1.4 M_sun, R = 5.8M: 1 + z = 1.2355; g = 0.03674/M = 1.6e12 m/s² (≈ 1.6e11 g).
- The Sun: z = 2.12e-6, a redshift of 0.636 km/s (measured: 633 m/s).
- Sirius B: R/M = 3,877 → z = 2.580e-4 → 77.3 km/s (measured 80.65 ± 0.77 km/s, with a
  slightly smaller radius).

**Rigidly rotating surface** with angular velocity Ω:

```
u^t = [ −(g_tt + 2Ω g_tφ + Ω² g_φφ) ]^(−1/2),   u^φ = Ω u^t
Frequency ratio for a photon with the code's normalization (camera energy 1):
g = 1 / ( u^t (E − Ω L) )          (same form as the disk in trace.frag)
```

**Light bending at the surface** (a = 0). A photon leaving the surface at angle α to the
radial direction (static frame) reaches infinity having swept an angle ψ around the star:

```
impact parameter      b = R sin α / √(1 − 2M/R)
exact                 ψ(α) = ∫_R^∞ dr / ( r² √( 1/b² − (1 − 2M/r)/r² ) )
Beloborodov (2002)    1 − cos α ≈ (1 − cos ψ)(1 − u),   u = 2M/R     (accurate to ~1 % for R > 4M)
Visible part of the surface: α = π/2 gives ψ_max.  Beloborodov: cos ψ_max = −u/(1 − u)
   → R = 5.8M: exact ψ_max = 120.67° (Beloborodov 121.76°); you see (1 − cos ψ_max)/2 = 75.5 % of the surface.
   → The whole surface is visible for R < 3.5206M (exact, from the integral; Beloborodov says 4M).
   (Numerics: substitute r = R/x, x ∈ (0, 1]; the α = π/2 end point is an integrable singularity.)
```

**Emission angle at the surface** in the surface's own frame (for limb darkening): build the
frame with `kerr::observerFrame(u_surface)` and take μ = (p·e_r) / (−p·u). For a = 0 and a
static surface this is μ = |vr| / (r² E).

**Limb darkening** (Eddington): I(μ)/I(1) = (2 + 3μ)/5. Use it for stars and white dwarfs.

### 2.5 How each object looks

- **Normal stars:** blackbody at T_eff with limb darkening, granulation (cells ~1,000 km on the
  Sun), sunspots and faculae, an optically thin corona above the surface, I ∝ exp(−(r − R)/H).
- **White dwarf:** smooth, hot, blue-white; optional magnetic spots.
- **Neutron star:** surface ~1e6 K (draw it with a visible color-temperature mapping and say
  so), two magnetic-pole hot spots (~1e7 K) on an axis tilted from the spin axis, the existing
  pulsar beams, and strong light bending (you see the back of the star). Magnetar option:
  surface flares.
- **Planet:** albedo map lit by the disk (reuse the existing planet shading), atmosphere rim.
- **Black hole:** as now.

### 2.6 Accretion onto a surface versus a horizon (the "no surface" evidence)

```
Black hole: energy released per unit rest mass  η = 1 − E_isco
            (5.72 % at a = 0, ~32 % at a = 0.998, from kerr::circularOrbit)
Star:       disk to the ISCO (or to the surface if R > r_isco) releases 1 − E_in;
            matter then settles on the surface at rest with E_surface = √(1 − 2M/R)
            (for a = 0), so the boundary layer and impact release E_in − E_surface.
            Neutron star R = 5.8M, a = 0: total 1 − √(1 − 2/5.8) = 19.1 %.
            White dwarf R/M = 9,800: ~1e-4.  Sun: ~2e-6.
Black hole: the same matter crossing the horizon releases nothing more: no impact flash.
```

Implement impacts: when a body or gas parcel reaches the surface, remove it, log
"hit the surface", and add a hot spot there whose brightness reflects the released energy and
decays over a few M. On a neutron star this is an X-ray-burst-like flash. **Type I X-ray
bursts** (thermonuclear, about 5 MeV per nucleon, so about 0.5 % of rest mass) can be an
occasional surface-wide flash. **Novae** do the same on a white dwarf.

**Magnetospheric truncation** (neutron stars, white dwarfs): the disk stops at the Alfvén radius

```
r_A = ( μ_B⁴ / (2 G M Ṁ²) )^(1/7),   μ_B = B R³ (dipole moment),   disk inner edge ≈ 0.5 r_A
Inside, gas flows along dipole field lines r = L sin²θ_m to the magnetic poles (accretion columns).
```

### 2.7 Neutron star and white dwarf structure (for the mass–radius relation and collapse)

```
TOV equations (G = c = 1):
  dm/dr = 4π r² ε
  dP/dr = −(ε + P)(m + 4π r³ P) / ( r (r − 2m) )
Polytrope: P = K ρ₀^Γ,  ε = ρ₀ + P/(Γ − 1)
Standard test (G = c = M_sun = 1, so lengths in 1.4766 km):
  K = 100, Γ = 2, ρ_c = 1.28e-3  →  M = 1.400 M_sun, circumferential R = 9.586 (14.15 km);
  the maximum mass for this equation of state is 1.637 M_sun.
  (Stop at P < 1e-14; these match the standard "TOV star" test of numerical relativity codes.)

White dwarf mass–radius (Nauenberg 1972, μ_e = 2):
  R ≈ 0.0112 R_sun √( (M/M_Ch)^(−2/3) − (M/M_Ch)^(2/3) ),   M_Ch = 1.44 M_sun
  → 0.6 M_sun gives R ≈ 8,660 km.
```

### 2.8 The compactness slider and collapse

- A log slider for R/M from about 5e5 (the Sun) down to 9/4 (the Buchdahl limit), mass fixed.
  Show the live readouts from 2.4 and 2.6, plus "orbits outside the surface: unchanged".
- Passing below the Buchdahl limit, or above the TOV maximum mass in neutron star mode, or
  above M_Ch in white dwarf mode, **collapses** the object into a black hole with the same M and J.
  If χ = J/M² ≥ 1, the collapse must shed angular momentum first: clamp to 0.998 and say so.
  Show a flash and an event message.
- The disk inner edge becomes max(ISCO, R) (or the Alfvén radius when magnetized). The photon
  ring appears only once R < 3M.
- Optional "true scale" versus "exaggerated" toggle: physics readouts stay real while the
  visuals can exaggerate compactness, because for the Sun all relativistic effects are ~1e-6.

### 2.9 Code changes

- New `src/central.h/.cpp`: `struct CentralObject { Kind kind; double massSolar; double radius;
  double spin; double surfaceTemperature; double surfaceOmega; vec3 magneticAxis;
  double magneticField; … }`, the catalogue, unit conversions and derived readouts (2.4, 2.6),
  `hasHorizon()` and `boundaryRadius()` (r+ for a black hole, R for a star).
- Settings and UI: rename the "Black hole" panel to "Central object". Add a kind selector with
  the presets from 2.2; sliders for mass (M_sun, log), radius (km or R/M; defaults from the
  relations in 2.7), spin (limited as in 2.3) and surface temperature; the compactness slider;
  the readouts.
- Replace every use of `kerr::horizonRadius` as a *boundary* with the central object's
  boundary (grep for `horizonRadius(`): ray termination, body removal, gas `rDead`, camera
  limits, the `classify` and `predict` fates ("hits the surface" versus "falls in"), the disk
  inner edge, and the UI tables. Keep the black-hole-only interior logic behind `hasHorizon()`.
- `trace.frag`: new uniforms for the central kind, R, T, Ω, magnetic axis, hot spots and impact
  spots. In the step loop, if `r1 ≤ R` for a star, interpolate the hit point, shade the surface
  (blackbody at g·T, limb darkening, texture in a frame rotating at Ω; for a ≠ 0 rotate with
  the Boyer–Lindquist azimuth φ_BL = φ − azimuthShift(r)), and stop the ray.
- The camera landing on a star: stop at the surface with a message giving the impact speed and
  the surface gravity in g.
- Scene scale for non-compact objects: default camera distance about 4R, log sliders up to
  ~1e9 M, escape radius above the camera distance, debris grid radii scaled to the object
  (the grid is in log r, so only its bounds change). Watch float32 precision in the shader at
  r ~ 1e6–1e9: either work in units of R on the GPU or accept straight-line rays there.
- **Orbiting compact objects:** add body kinds WhiteDwarf, NeutronStar (non-pulsing) and
  MiniBlackHole. A mini black hole cannot be ray traced exactly inside a Kerr background (no
  exact two-black-hole spacetime exists). Approximate it as a thin gravitational lens: a ray
  passing with impact parameter b is bent by 4m/b toward it, applied as a momentum kick to the
  photon (recompute E, L, Q from the new direction, as gravity kicks do for bodies), and is
  captured if b < 3√3 m. Document the approximation.
- **Tidal disruption needs the right scales:** r_t = R* (M/m*)^(1/3). A star is **swallowed
  whole** if r_t < r+ (more precisely, the marginally bound radius).

```
Hills mass  M_H = ( R* c² / (2G) )^(3/2) / m*^(1/2)
  Sun-like star: M_H ≈ 1.1e8 M_sun.  Typical white dwarf: ≈ 2e5 M_sun.
In the simulator's units: R*/M = R*[km] / (1.4766 km × M/M_sun).
  Sgr A*: the Sun's radius is 0.11 M → r_t ≈ 18 M (disrupted, visible).
  M87*:   r_t ≈ 0.14 M < r+ → swallowed whole.
```

---

## 3. The ten other features

### 3.1 Penrose process

A body enters the ergosphere and splits. One piece falls in with negative energy; the other
escapes with more energy than the original had. The energy comes out of the hole's rotation.

```
ZAMO-frame quantities: α = √(ΣΔ/A),  ω = 2ar/A,  ϖ = √(A/Σ) sinθ
For ZAMO-frame velocity v (γ = 1/√(1 − v²)):  E = γ(α + ω ϖ v_φ),  L = γ ϖ v_φ   (per unit mass)
E < 0 needs v_φ < −α/(ω ϖ), possible with |v| < 1 only where α < ω ϖ: inside the ergosphere.

Split (local 4-momentum conservation, in ZAMO orthonormal components):
  m0 γ0 (1, v0) = m1 γ1 (1, v1) + m2 γ2 (1, v2),   m1 + m2 ≤ m0
  → m0 E0 = m1 E1 + m2 E2,   m0 L0 = m1 L1 + m2 L2
Choose fragment 1 counter-rotating with E1 < 0; solve P2 = P0 − P1; require m2 = √(−P2·P2) real
and P2 future-directed.
Energy gained:  m2 E2 − m0 E0 = −m1 E1 > 0
Hole after fragment 1 crosses the horizon:  δM = m1 E1 < 0,  δJ = m1 L1
Area theorem check: δM − Ω_H δJ = m1 (E1 − Ω_H L1) ≥ 0  with Ω_H = a / (2 r+)
   (this is exactly the P(r+) ≥ 0 condition the code already uses for horizon crossings)
Maximum efficiency (a → 1, split at the horizon): (√2 − 1)/2 ≈ 20.7 %
Irreducible mass: M_irr = √(A/16π) = √(r+/2) (M = 1);  M² = M_irr² + J²/(4 M_irr²)  (Christodoulou)
Extractable rotational energy: M − M_irr  (29 % at a = 1)
```

UI: a "Penrose process" preset (equatorial, a ≥ 0.9), a panel showing E0, E1, E2, efficiency,
and the hole's M, J, a, M_irr and area before and after. Optionally apply δM and δJ to the hole
(spin down) through the existing respin path.

### 3.2 Gravitational-wave inspiral with a chirp sound

A small body (a mini black hole or a neutron star of mass μ) loses energy and angular momentum
to gravitational waves, spirals in and plunges.

```
Circular, equatorial (adiabatic sequence of circular orbits), M = 1, q = a:
  v = (Ω)^(1/3),   Ω = 1 / (r^(3/2) + q)
  Energy flux (Tagoshi, Shibata, Tanaka & Sasaki 1996, test-mass limit):
    F = (32/5) μ² v^10 [ 1 − (1247/336) v² + (4π − (11/4) q) v³
                         + (−44711/9072 + (33/16) q²) v⁴ + (−(8191/672) π − (59/16) q) v⁵ + … ]
  Orbital energy per unit μ:  E(x) = (1 − 2x² + q x³) / √(1 − 3x² + 2q x³),   x = √(1/r)
  Radius evolution:  dr/dt = −(F/μ) / (dE/dr)      → diverges at the ISCO: switch to a geodesic
                                                     plunge (Ori & Thorne 2000 transition).
  Orbital phase:  dφ/dt = Ω;  GW phase Φ = 2φ;  GW frequency f_GW = Ω/π  (in 1/M; Hz = /T_M)
  Strain at distance D:  h₊ = −(4μ/D)(Ω)^(2/3) (1 + cos² ι)/2 cos Φ,   h× = −(4μ/D)(Ω)^(2/3) cos ι sin Φ

Eccentric orbits (Peters & Mathews 1963; a = semi-major axis, e = eccentricity, M = total mass):
  ⟨dE/dt⟩ = −(32/5) μ² M³ / a⁵ · (1 + 73/24 e² + 37/96 e⁴) / (1 − e²)^(7/2)
  ⟨dL/dt⟩ = −(32/5) μ² M^(5/2) / a^(7/2) · (1 + 7/8 e²) / (1 − e²)²
  ⟨da/dt⟩ = −(64/5) μ M² / a³ · (1 + 73/24 e² + 37/96 e⁴) / (1 − e²)^(7/2)
  ⟨de/dt⟩ = −(304/15) e μ M² / a⁴ · (1 + 121/304 e²) / (1 − e²)^(5/2)
  Circular coalescence time:  t_c = (5/256) a⁴ / (μ M²);   a(t)⁴ = a₀⁴ − (256/5) μ M² t
```

Implementation for generic orbits ("osculating geodesics"): each body with mass μ > 0 gets
radiation reaction. Get p and e from the turning points (`trajectory::classify`) and the
inclination ι = atan(√Q / |L|). Each step, change E and L by the fluxes, keep ι (Q = L² tan² ι),
then put the velocities back on vr² = R, vθ² = Θ, keeping their signs. Add a
**"radiation-reaction boost"** slider (real EMRIs take ~10⁵ orbits) and show the true time to
merger. The plunge adds μ to the hole's mass. Around a star, the inspiral ends in an impact.

**Sound:** a 10 M_sun black hole has, at a = 0, f_GW at the ISCO = 1/(π 6^(3/2) T_M) ≈ **440 Hz**
(concert A; make this a test). For physical masses play the true f_GW. For supermassive holes
(mHz, LISA band) scale the pitch and show the true frequency. Synthesize
amplitude ∝ Ω^(2/3) and the phase by integrating 2Ω.

### 3.3 "What the telescope sees" (Event Horizon Telescope mode)

```
Radiative transfer along each backward ray (camera energy normalized to 1, g = 1/(−p·u)):
  fluid-frame path length per affine step:   dℓ_f = dλ / g,   with dλ = Σ dτ_Mino
                                             (in units of M; multiply by L_M in cm for CGS j and α)
  invariant intensity 𝓘 = I_ν / ν³:          d𝓘/dλ = (1/g) [ j_ν(ν_f)/ν_f³ − α_ν(ν_f) 𝓘 ],   ν_f = ν_obs / g
  front-to-back:  I_obs = ν_obs³ Σ_k exp(−τ_k) (j_k / ν_f,k³) Δλ_k / g_k,   τ_{k+1} = τ_k + α_k Δλ_k / g_k

Thermal synchrotron emissivity (Leung, Gammie & Noble 2011), CGS:
  ν_c = e B / (2π m_e c) = 2.80e6 Hz × B[G],   θ_e = k T_e / (m_e c²)
  ν_s = (2/9) ν_c θ_e² sin θ_B,   X = ν / ν_s
  j_ν = n_e (√2 π e² ν_s) / (3 K₂(1/θ_e) c) · (X^(1/2) + 2^(11/12) X^(1/6))² exp(−X^(1/3))
        (K₂(1/θ_e) ≈ 2 θ_e² for θ_e ≳ 1;  θ_B is the angle between B and the photon in the fluid frame)
  Absorption (Kirchhoff):  α_ν = j_ν / B_ν(T_e),   B_ν = (2hν³/c²) / (exp(hν/kT_e) − 1)

Hot "puffy" flow (radiatively inefficient, Broderick & Loeb style; tune the constants):
  n_e = n₀ (r/M)^(−1.1) exp(−z² / (2ρ²)),   z = r cos θ,  ρ = r sin θ
  T_e = T₀ (r/M)^(−0.84)
  B   = √( 8π n_e m_p c² / (6 β (r/M)) ),   β ≈ 10
  Gas motion: sub-Keplerian rotation outside the ISCO, geodesic infall inside.
  Calibrate n₀ so the total flux at 230 GHz is ≈ 0.5 Jy (M87*) or ≈ 2.4 Jy (Sgr A*).

Observer at infinity (Bardeen 1973), image-plane coordinates (α, β) in units of M, inclination i:
  L/E = −α sin i,   Q/E² = β² + (α² − a²) cos² i,   sign of p_θ from the sign of β
Angular scale: θ_g = G M / (c² D) → M87* (6.5e9 M_sun, 16.8 Mpc): 3.82 µas;  Sgr A* (4.3e6, 8.28 kpc): 5.13 µas
Flux:  F_ν [Jy] = 1e23 Σ_pixels I_ν Δα Δβ θ_g²;   brightness temperature T_b = c² I_ν / (2 k ν²)
Blur: Gaussian, FWHM 20 µas (σ = FWHM / 2.3548).  Color map: "afmhot".
Shadow edge (critical curve) for any a, parametrized by the photon-orbit radius r:
  ξ(r) = L/E = [r²(3 − r) − a²(r + 1)] / [a (r − 1)],   η(r) = Q/E² = r³ [4a² − r(r − 3)²] / [a² (r − 1)²]
  α = −ξ / sin i,   β = ±√(η + a² cos² i − ξ² cot² i)
a = 0 shadow diameter: 2√27 M = 10.39 M → M87*: 39.7 µas (EHT ring 42 ± 3 µas);  Sgr A*: 53.3 µas (51.8 ± 2.3)
```

Show the simulated image next to the real EHT images (link to them; do not embed copyrighted
images) with the same angular scale. Polarization (EVPA ticks) can come later.

### 3.4 A disk that grows from tidal debris

```
Most bound debris: specific energy spread Δε ≈ G M R* / r_t²  →  a_min = r_t² / (2 R*)
                   t_min = 2π √(a_min³ / (G M))
Fallback rate (flat dM/dε):  Ṁ(t) = (m* / (3 t_min)) (t / t_min)^(−5/3),  t > t_min
Circularization radius (Newtonian, L conserved): r_c = 2 r_p   (in Kerr: the circular orbit with the star's L)

1-D disk, surface density Σ(r, t):
  ∂Σ/∂t = (3/r) ∂/∂r [ r^(1/2) ∂/∂r ( ν Σ r^(1/2) ) ] + S(r, t)
  α-viscosity: ν = α c_s H = α c_s² / Ω  (α ≈ 0.1);  inner boundary Σ = 0 at r_in (ISCO), free outer edge
  explicit step limit: Δt ≲ Δr² / (12 ν)  (conservative; or use an implicit scheme)
  Newtonian local dissipation: D(r) = (9/8) ν Σ Ω²,   σ T⁴ = D / 2 per face
  Steady state check: ν Σ = (Ṁ / 3π)(1 − √(r_in / r))

Relativistic steady-state flux (Page & Thorne 1974), M = 1, x = √r, x₀ = √r_isco:
  x₁ = 2 cos( (1/3) arccos a − π/3 ),  x₂ = 2 cos( (1/3) arccos a + π/3 ),  x₃ = −2 cos( (1/3) arccos a )
  F(r) = 3Ṁ / (8π x⁴ (x³ − 3x + 2a)) · [ x − x₀ − (3/2) a ln(x/x₀)
          − 3(x₁ − a)² / (x₁ (x₁ − x₂)(x₁ − x₃)) · ln((x − x₁)/(x₀ − x₁))
          − 3(x₂ − a)² / (x₂ (x₂ − x₁)(x₂ − x₃)) · ln((x − x₂)/(x₀ − x₂))
          − 3(x₃ − a)² / (x₃ (x₃ − x₁)(x₃ − x₂)) · ln((x − x₃)/(x₀ − x₃)) ]
  (at a = 0, x₂ = 0 and its term vanishes: drop it for |a| < 1e-6)
  T(r) = (F / σ)^(1/4)
```

Connect it to the gas: when a parcel's orbit is nearly circular and near the disk plane
(e < 0.1, |cos θ| small) or has thermalized in a shock, remove it from the particle system and
add its mass to the disk ring at its radius. The renderer's disk then uses T(r) and opacity from
Σ (upload a 1-D texture) instead of the fixed profile. The disk lights up over time like a real
tidal-disruption flare.

### 3.5 Naked singularity (spin above 1)

```
For a > 1: Δ = r² − 2r + a² = (r − 1)² + b² > 0 for all r, with b = √(a² − 1). No horizons.
Ring singularity at r = 0, θ = π/2 (Σ = 0). Rays through the ring's disk (r = 0, θ ≠ π/2) continue
to r < 0: another asymptotic region where gravity repels. Mino-time R(r) and Θ(θ) are regular at
r = 0, so the integration simply continues through it.
Kerr–Schild shifts in closed form:
  ∫ (r² + a²)/Δ dr − r = ln Δ + (2/b) arctan((r − 1)/b)      (derivative: 2r/Δ)
  ∫ a/Δ dr             = (a/b) arctan((r − 1)/b)              (derivative: a/Δ)
Equatorial photon orbits: r (r − 3)² = 4a²  (solve numerically; the closed forms using
acos(±a) fail for a > 1).  ISCO: solve dE/dr = 0 for circular orbits numerically.
Shadow: the critical curve exists only partly → an "arc" (Hioki & Maeda 2009).
```

Code: branch everything that uses √(1 − a²) on `hasHorizon()`. The Cartesian chart maps r < 0
onto the same region as r > 0, so carry a "sheet" flag. Draw the r → −∞ sky differently (for
example tinted or mirrored), and draw rays that hit Σ < ε as the singularity. Add a warning in
the UI: cosmic censorship, and ergoregion instability without a horizon.

### 3.6 Timeline scrubbing

- Make the simulation deterministic: fixed simulation ticks (say 0.25 M), seeded random number
  generators (BodySystem already has one), deterministic GPU updates. Integer fixed-point
  `atomicAdd` grids are already order-independent.
- Checkpoints every ~25 M for the last ~2,000 M: the CPU state (bodies, departed bodies, events,
  time, RNG state, camera) and the GPU particle buffer, copied into a ring of buffers sized by
  the active particle count, within a memory budget (~512 MB).
- The light-delay history (`History` in `src/bodies.h`: `kColumns` = 4096 samples, `kStep` = 0.5 M,
  `newest`) is indexed by absolute time, and a deterministic re-run writes the same samples. After
  restoring to time T, set `newest = floor(T / kStep)`. Careful: samples from the abandoned future
  sit in the ring slots of the oldest times, so either store the ring in each checkpoint (it is
  small) or treat that range as missing until it is overwritten.
- Scrubbing to T means restoring the latest checkpoint before T and re-simulating forward
  without rendering. For reverse play, step back in chunks.
- UI: a timeline bar (buffered range, current time, play / pause / reverse / speed). Keys: `[`
  and `]` step, Shift+Space plays backwards.
- Test: run, checkpoint, run on, restore, re-run → identical states (bitwise for the CPU, read
  back for the GPU).

### 3.7 Sound

- Use miniaudio (single header, public domain or MIT-0; fetch it with CMake FetchContent). Run
  the audio callback on its own thread; pass parameters through atomics or a lock-free ring;
  smooth them (one-pole filters) to avoid clicks.
- **Every sound is carried by light:** its frequency and rate are multiplied by g, the same
  frequency ratio the ray tracer uses. Get g and the delay per body from a CPU ray toward it
  (`trajectory::traceRay`) a few times a second.
- **Pulsar:** a click each time the beam sweeps the line of sight at emission time
  (t_emit = t_now − delay); click rate = spin rate × g; loudness ∝ g³.
- **Bodies:** a soft tone at f₀ × g, louder when close, panned by direction in the camera frame.
- **Near the hole:** a low rumble with amplitude ∝ the tidal acceleration at the camera, ~2M/r³.
- **Gravitational waves:** the chirp from 3.2.
- Volume and mute controls; keep it quiet by default.

### 3.8 Cinematic camera keyframes

```
Keyframe: Cartesian position (Kerr–Schild), orientation quaternion, field of view, optional hold.
Position: centripetal Catmull–Rom (α = 0.5); knots t_{i+1} = t_i + |P_{i+1} − P_i|^α; Barry–Goldman evaluation.
Constant speed: reparametrize by arc length (numeric integration + lookup table), then ease
  with smoothstep s(u) = u² (3 − 2u).
Orientation: squad for C¹ smoothness:
  squad(q_i, q_{i+1}, s_i, s_{i+1}; h) = slerp( slerp(q_i, q_{i+1}, h), slerp(s_i, s_{i+1}, h), 2h(1 − h) )
  s_i = q_i exp( −( log(q_i⁻¹ q_{i+1}) + log(q_i⁻¹ q_{i−1}) ) / 4 )
```

Two modes: **cinema** (the camera is at rest relative to the local ZAMO at every point, so no
aberration) and **physical** (the 4-velocity comes from the path's velocity, giving real
aberration; speed must stay below c and the path must stay outside the horizon). UI: add a
keyframe (K), list, reorder, set durations, preview, record.

### 3.9 Save and load scenes

- JSON with a version field (nlohmann/json via FetchContent, or a small writer). Store: the
  central object, settings, camera (including its geodesic), every body (constants, state, mass,
  radius, kind, temperature, shape F and Ḟ, maxStretch, inside), departed bodies, events, time,
  RNG seeds. Store debris in a binary sidecar (`.bhscene.bin`, the float32 particle array).
- Save with Ctrl+S, open with Ctrl+O, drag and drop a file onto the window (`glfwSetDropCallback`),
  and a `--scene-file PATH` flag.
- Test: save, load, save again → identical files and identical simulation afterwards.

### 3.10 Speed

- **Dynamic resolution:** scale_{n+1} = clamp( scale_n · (target_ms / measured_ms)^0.5, 0.4, 1.0 ),
  smoothed.
- **FSR 1** (AMD FidelityFX Super Resolution 1, MIT licence; EASU upscaling + RCAS
  sharpening): run it after tone mapping.
- **Static-camera cache:** when the camera and spin are still, the photon paths are identical
  frame to frame. Cache per pixel the ray's end (sky direction, or disk hits with r, φ, g, delay
  and transmittance) and re-shade the moving disk pattern and sky without re-tracing. Re-trace
  fully only pixels near moving content (bodies, gas, paths, jets), using a mask from the
  previous frame dilated by the lensing spread, and always the photon-ring band.
- GPU timer queries per pass (trace, bloom, particles) shown in the stats line.

---

## 4. Tests to add (closed-form or published values)

1. Static surface redshift at R = 5.8M: a radial photon arrives with g = √(1 − 2/5.8) (to 1e-9).
2. Surface light bending at R = 5.8M: the traced ψ(α) matches the exact integral (to 1e-6 rad)
   for several α; ψ_max = 120.67°; the whole surface is visible below R = 3.5206M.
3. Mercury: the exact Schwarzschild precession formula (already in the tests) at
   p = a(1 − e²) with a = 5.7909e10 m, e = 0.20563 gives 42.98″ per century (415.2 orbits).
4. Sun: limb deflection 4M/b = 1.751″; surface redshift 0.636 km/s. Sirius B: 77.3 km/s.
5. TOV (K = 100, Γ = 2, ρ_c = 1.28e-3): M = 1.4002 M_sun, R = 9.586; maximum mass 1.637 M_sun
   (near ρ_c = 3.2e-3). Nauenberg: 0.6 M_sun → 8,658 km.
6. Accretion efficiency: 1 − √(8/9) = 5.72 % (black hole, a = 0); 19.1 % (neutron star R = 5.8M).
7. Hills mass for a Sun-like star: 1.14e8 M_sun.
8. Penrose: conservation of E and L at the split; E1 − Ω_H L1 ≥ 0 for every infalling fragment;
   Christodoulou's identity M² = M_irr² + J²/(4 M_irr²).
9. Gravitational waves: the Peters a(t) law against the integrated inspiral at large r; the leading
   PN flux term; f_GW at the ISCO for 10 M_sun = 439.7 Hz.
10. EHT: the rendered shadow edge matches the critical curve (a = 0.9, i = 17°); diameter
    2√27 M at a = 0; the M87* shadow ≈ 39.7 µas.
11. Disk: the steady-state νΣ profile; Page–Thorne F(r_isco) = 0 and F → 3Ṁ/(8π r³) as r → ∞.
12. Naked singularity: Δ > 0; derivatives of the closed-form shifts; equatorial photon orbits
    satisfy r(r − 3)² = 4a²; a ray through the ring's disk keeps E, L and Q.
13. Keyframes: the spline passes through every keyframe; arc-length speed varies < 1 %.
14. Save and load round trip; timeline determinism.

---

## 5. Order of work and how to report

Work in phases. Finish each one completely (code, tests, screenshots checked by eye, GPU time
measured, README updated) before starting the next, and give me a short report after each:
what changed, test results, timings, screenshots, and what is approximate.

1. Central object selector, surfaces, physical units, compactness slider and collapse (section 2).
2. Penrose process (3.1).
3. Gravitational-wave inspiral and the chirp (3.2), with sound (3.7).
4. Timeline scrubbing (3.6) and save/load (3.9).
5. Event Horizon Telescope mode (3.3).
6. Disk from debris (3.4).
7. Naked singularity (3.5).
8. Camera keyframes (3.8).
9. Speed (3.10).

Do not break existing features: every current preset, scene and flag must still work, and all
existing tests must still pass. Ask me before big refactors, before deleting anything, and
before pushing to GitHub.
