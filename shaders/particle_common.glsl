// Shared by the debris compute shaders: the particle layout and the geodesic
// equations for a massive particle (mu = 1), the GPU twin of kerr::derivatives,
// plus conversions to and from the local ZAMO frame.
//
// The GPU works in render units of `scale` M (src/central.h, Object::sceneScale): 1 for
// black holes and neutron stars, R/6 around bigger objects, so numbers stay of order
// 1-1000 in single precision. The metric is then Kerr with mass m = 1 / scale and spin
// a = chi m; every formula carries m explicitly (m = 1 is the original M = 1 form).

struct Particle {
    vec4 state;      // r, theta, phi, and 1 = alive, 0.5 = gone but its light still on the way, 0 = dead
    vec4 velocity;   // vr, vtheta (Mino time), temperature (K), weight
    vec4 constants;  // E, L, Q, entropy T / rho^(gamma - 1) (0 until the gas around it is known)
    // The parcel as seen now: where it was when the light now reaching the camera left
    // it (light-travel delay). r, theta, phi, and that time. Seen only for times after
    // it was created (-1e30 if it has always been there) and before it was lost.
    vec4 seen;
    vec4 seenVelocity;  // vr, vtheta, time created, time lost (1e30 while it lasts)
};
layout(std430, binding = 1) buffer Particles {
    Particle particles[];
};

const float PI = 3.14159265358979;

// (P + vr) / Delta without the 0/0 at the horizon (see kerr::inflowFactor).
float inflowFactor(float P, float vr, float delta, float K) {
    if (P * vr < 0.0) return K / (P - vr);
    return (P + vr) / (abs(delta) > 1e-6 ? delta : (delta < 0.0 ? -1e-6 : 1e-6));
}

// Mino-time derivatives in Kerr-Schild coordinates: dt/dtau, dphi/dtau, d2r/dtau2,
// d2theta/dtau2, for radial velocity vr = dr/dtau, around mass m with spin a.
void massiveRhs(float m, float a, float E, float L, float Q, float r, float th, float vr,
                out float tDot, out float phiDot, out float accR, out float accTh) {
    float s = sin(th), c = cos(th);
    float s2 = max(s * s, 1e-8);
    float s3 = s2 * max(abs(s), 1e-4) * (s < 0.0 ? -1.0 : 1.0);
    float a2 = a * a;
    float delta = r * r - 2.0 * m * r + a2;
    float rra = r * r + a2;
    float P = E * rra - a * L;
    float x = L - a * E;
    float K = r * r + x * x + Q;
    float F = inflowFactor(P, vr, delta, K);
    accR = 2.0 * E * r * P - (r - m) * K - delta * r;
    accTh = s * c * a2 * (1.0 - E * E) + L * L * c / s3;
    phiDot = a * F - a * E + L / s2;
    tDot = rra * F - a * (a * E * s * s - L) - vr;
}

// Re-projects the velocities onto vr^2 = R(r) and vtheta^2 = Theta(theta) away from
// turning points (see kerr::enforceConstraints).
void enforceMassiveConstraints(float m, float a, float E, float L, float Q, float r, float th, inout float vr,
                               inout float vth) {
    float a2 = a * a;
    float delta = r * r - 2.0 * m * r + a2;
    float P = E * (r * r + a2) - a * L;
    float x = L - a * E;
    float K = r * r + x * x + Q;
    float R = P * P - delta * K;
    if (R > 0.05 * (P * P + abs(delta) * K)) vr = (vr < 0.0 ? -1.0 : 1.0) * sqrt(R);

    float s = sin(th), c = cos(th);
    float cot2 = c * c / max(s * s, 1e-8);
    float polarMass = a2 * (1.0 - E * E) * c * c;
    float Theta = Q - polarMass - L * L * cot2;
    if (Theta > 0.05 * (abs(Q) + abs(polarMass) + L * L * cot2)) vth = (vth < 0.0 ? -1.0 : 1.0) * sqrt(Theta);
}

// The outer horizon r+ = m + sqrt(m^2 - a^2).
float horizonRadius(float m, float a) {
    return m + sqrt(max(m * m - a * a, 0.0));
}

// --- The local ZAMO frame (outside the horizon only) ------------------------------

struct Zamo {
    float sigma, delta, bigA, omega;
    float alpha;  // lapse: the ZAMO's clock rate
    float varpi;  // sqrt(A / Sigma) sin(theta): circumference radius
};

Zamo zamoAt(float m, float a, float r, float th) {
    float s = max(abs(sin(th)), 1e-6), c = cos(th);
    Zamo z;
    z.sigma = r * r + a * a * c * c;
    z.delta = r * r - 2.0 * m * r + a * a;
    float rra = r * r + a * a;
    z.bigA = rra * rra - a * a * z.delta * s * s;
    z.omega = 2.0 * m * a * r / z.bigA;
    z.alpha = sqrt(max(z.sigma * z.delta / z.bigA, 1e-12));
    z.varpi = sqrt(z.bigA / z.sigma) * s;
    return z;
}

// The particle's 3-velocity in the ZAMO frame, (r, theta, phi) components (kerr::zamoVelocity).
vec3 zamoVelocity(Zamo z, float E, float L, float vr, float vth) {
    float gamma = (E - z.omega * L) / z.alpha;
    return vec3(vr / sqrt(z.sigma * z.delta), vth / sqrt(z.sigma), L / z.varpi) / gamma;
}

// Constants of motion and Mino velocities for ZAMO-frame velocity v (kerr::launchFromZamo).
void launchMassive(float a, Zamo z, float th, vec3 v, out float E, out float L, out float Q, out float vr,
                   out float vth) {
    float s = max(abs(sin(th)), 1e-6), c = cos(th);
    float v2 = dot(v, v);
    if (v2 > 0.998) v *= sqrt(0.998 / v2);
    float gamma = inversesqrt(1.0 - dot(v, v));
    L = gamma * v.z * z.varpi;
    E = gamma * z.alpha + z.omega * L;
    vr = gamma * v.x * sqrt(z.sigma * z.delta);
    vth = gamma * v.y * sqrt(z.sigma);
    Q = vth * vth + c * c * (a * a * (1.0 - E * E) + L * L / (s * s));
}

// Kerr-Schild Cartesian position (kerr::toCartesian).
vec3 cartesian(float a, float r, float th, float ph) {
    float s = sin(th), cp = cos(ph), sp = sin(ph);
    return vec3((r * cp - a * sp) * s, (r * sp + a * cp) * s, r * cos(th));
}
