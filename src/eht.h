#pragma once

// What the Event Horizon Telescope sees: the black hole imaged at 230 GHz from Earth.
//
// Each pixel of the image plane of a distant observer (Bardeen 1973: impact parameters
// alpha, beta in M) is traced back into the hole's spacetime, and the light of a hot,
// "puffy" accretion flow (a radiatively inefficient flow, as around M87* and Sgr A*) is
// gathered along it by radiative transfer: thermal synchrotron emission (Leung, Gammie
// and Noble 2011) and its absorption (Kirchhoff), each seen at the frequency the gas emits
// at, nu_f = nu / g, Doppler and gravitationally shifted. The total flux sets the gas's
// density (n0), calibrated to what is observed (0.5 Jy for M87*, 2.4 Jy for Sgr A*). The
// image is then blurred to the array's resolution (20 microarcseconds).
//
// Units: G = c = M = 1 for the geometry; the plasma in cgs (cm^-3, K, gauss); intensities
// in erg s^-1 cm^-2 Hz^-1 sr^-1; fluxes in janskys.

#include <atomic>
#include <string>
#include <vector>

namespace eht {

// --- Physical constants (cgs) -------------------------------------------------------------

constexpr double kElectronCharge = 4.80320425e-10;  // esu
constexpr double kElectronMass = 9.1093837015e-28;  // g
constexpr double kProtonMass = 1.67262192369e-24;   // g
constexpr double kLight = 2.99792458e10;            // cm/s
constexpr double kBoltzmann = 1.380649e-16;         // erg/K
constexpr double kPlanck = 6.62607015e-27;          // erg s
constexpr double kSolarLengthCm = 1.476625e5;       // G M_sun / c^2
constexpr double kParsecCm = 3.0856775814913673e18;
constexpr double kJansky = 1e-23;                   // erg s^-1 cm^-2 Hz^-1
constexpr double kMicroarcsecond = 3.14159265358979323846 / 180.0 / 3600.0 / 1e6;  // rad

// --- The targets ---------------------------------------------------------------------------

struct Target {
    const char* name;
    double massSolar;
    double distancePc;
    double inclinationDeg;   // of the spin axis to the line of sight
    double flux;             // Jy at 230 GHz
    double ringUas, ringErrorUas;  // the ring's diameter as the EHT measured it
    const char* reference;
    const char* url;
};
// M87* (EHT Collaboration 2019, Paper I and VI) and Sgr A* (EHT Collaboration 2022, Paper I).
const std::vector<Target>& targets();

// theta_g = G M / (c^2 D), the angle M subtends, in microarcseconds.
double angularScale(double massSolar, double distancePc);

// --- The shadow -------------------------------------------------------------------------------

// The critical curve, the edge of the shadow, for spin a seen at inclination i (radians):
// for each spherical photon orbit radius r between the prograde and retrograde ones,
//   xi = L/E = [r^2 (3 - r) - a^2 (r + 1)] / [a (r - 1)],
//   eta = Q/E^2 = r^3 [4 a^2 - r (r - 3)^2] / [a^2 (r - 1)^2],
//   alpha = -xi / sin i,  beta = +-sqrt(eta + a^2 cos^2 i - xi^2 cot^2 i),
// kept where the root is real. For a = 0, a circle of radius sqrt(27). Points in M,
// going once round.
struct Point {
    double alpha, beta;
};
std::vector<Point> criticalCurve(double a, double inclination, int samples);
// Its radius seen from the image's origin in the direction psi (from the alpha axis).
double criticalRadius(const std::vector<Point>& curve, double psi);

// --- Synchrotron radiation ----------------------------------------------------------------

// Modified Bessel functions of the second kind (Abramowitz and Stegun 9.8.5-9.8.8).
double besselK0(double x);
double besselK1(double x);
double besselK2(double x);
// Thermal synchrotron emissivity (Leung, Gammie and Noble 2011, their eq. 72) at frequency
// nu (Hz) from electrons of density ne (cm^-3) at thetaE = k T_e / (m_e c^2) in a field B
// (G) at angle thetaB to the line of sight:
//   nu_s = (2/9) nu_c thetaE^2 sin(thetaB), nu_c = e B / (2 pi m_e c), X = nu / nu_s,
//   j = ne sqrt(2) pi e^2 nu_s / (3 K2(1/thetaE) c) (X^1/2 + 2^(11/12) X^1/6)^2 exp(-X^1/3).
// erg s^-1 cm^-3 Hz^-1 sr^-1.
double synchrotronEmissivity(double nu, double ne, double thetaE, double B, double sinThetaB);
// Planck's B_nu(T), erg s^-1 cm^-2 Hz^-1 sr^-1.
double planck(double nu, double T);
// Brightness temperature c^2 I / (2 k nu^2) of an intensity I at nu.
double brightnessTemperature(double intensity, double nu);

// --- The flow --------------------------------------------------------------------------------

// A hot, thick flow (after Broderick and Loeb 2006):
//   n_e = n0 r^-1.1 exp(-z^2 / 2 (h rho)^2),  z = r cos(theta), rho = r sin(theta), h = H/R
//   T_e = T0 r^-0.84
//   B = sqrt(8 pi n_e m_p c^2 / (6 beta r))  (a plasma beta of ~10 against the ions' virial pressure)
// rotating at a fraction of the Keplerian rate outside the ISCO and falling in along the
// ISCO's plunging geodesic inside it, where it thins out as it speeds up: the mass flowing
// in, ~ r^2 n u^r, stays the same (continuity).
struct Flow {
    double n0 = 1e5;          // cm^-3 at r = 1 M (calibrated to the flux)
    double T0 = 4e11;         // K at r = 1 M
    double beta = 10.0;
    double keplerian = 0.9;   // the rotation, as a fraction of Keplerian
    double height = 1.0;      // H/R, the flow's thickness
    double sinThetaB = 0.8165;  // a tangled field: sqrt(<sin^2>) over random directions, sqrt(2/3)
};
struct Plasma {
    double ne, Te, B;
};
Plasma plasmaAt(const Flow& flow, double r, double theta, double a = 0.0);
// The gas's 4-velocity (Boyer-Lindquist, contravariant: u^t, u^r, u^theta, u^phi).
void flowVelocity(const Flow& flow, double a, double r, double theta, double u[4]);

// --- The image ---------------------------------------------------------------------------------

struct View {
    double a = 0.9;
    double inclination = 0.3;  // radians
    double halfWidth = 16.0;   // the field's half-width, in M
    int size = 128;            // pixels across
    double frequency = 230e9;  // Hz
    double massSolar = 6.5e9;  // sets the plasma's path lengths (L_M in cm)
};

struct Ray {
    double intensity = 0.0;  // I_nu at the observer
    bool captured = false;   // fell into the hole
    double minRadius = 0.0;
};
// Traces the ray arriving at (alpha, beta) back from r = 1000 M, gathering the flow's light
// (withLight false: only follows it, to find the shadow).
Ray traceRay(const View& camera, const Flow& flow, double alpha, double beta, bool withLight = true);

struct Image {
    int size = 0;
    double halfWidth = 0.0;           // M
    std::vector<double> intensity;    // I_nu, row by row from the top (beta = +halfWidth)
    std::vector<unsigned char> shadow;  // 1 where the ray fell in
    double pixel() const { return 2.0 * halfWidth / size; }  // M per pixel
};
// Renders on `threads` threads; stops early (returning what it has) if cancel turns true.
Image render(const View& camera, const Flow& flow, int threads, const std::atomic<bool>* cancel = nullptr);

// Total flux (Jy) of an image of a source whose M subtends thetaGUas microarcseconds.
double flux(const Image& image, double thetaGUas);
// A Gaussian blur of the given full width at half maximum (M).
Image blur(const Image& image, double fwhm);
// The ring's diameter (M): twice the mean, over 72 directions from (alpha, beta) =
// center, of the radius where the brightness peaks between 0.4 and 2 times `radius`.
double ringDiameter(const Image& image, Point center, double radius);
// The middle of a closed curve (the mean of its points).
Point curveCenter(const std::vector<Point>& curve);
// The density n0 that gives the observed flux (secant iteration in log n0 on small images).
double calibrate(const View& camera, Flow flow, double targetJy, double thetaGUas, int threads,
                 const std::atomic<bool>* cancel = nullptr);

// matplotlib's "afmhot" color map, as the EHT's images use: x in [0, 1] to RGB in [0, 1].
void afmhot(double x, float rgb[3]);

}  // namespace eht
