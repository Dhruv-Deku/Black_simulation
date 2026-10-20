#include "film.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace film {

namespace {

// The first shot's stars, tracked in every frame of the reference (the two brightest blobs,
// the nearer one the larger), smoothed by local quadratic fits and sampled every 0.1 s from
// 0 to 5.8 s. kPhase: the orbit's phase psi in radians (it turns 6 times a second at first,
// 30 times as they touch); kSemiMajor: the projected orbit's tall semi-axis, reference
// pixels; kGlare: the radius of a star's white glare where both are equally far, pixels.
constexpr double kStep = 0.1;
constexpr std::array<double, 59> kPhase = {
    3.095,  3.510,  4.109,  4.705,  5.219,  5.765,  6.318,  7.002,  7.553,  8.136,  8.820,  9.442,
    9.964,  10.531, 11.232, 11.932, 12.441, 13.061, 13.721, 14.306, 14.926, 15.570, 16.135, 16.744,
    17.540, 18.390, 19.070, 19.916, 20.920, 21.743, 22.650, 23.836, 24.881, 25.986, 27.245, 28.423,
    29.487, 30.842, 32.201, 33.686, 34.950, 36.468, 37.938, 39.529, 41.145, 42.688, 44.469, 46.571,
    48.051, 49.813, 52.102, 54.248, 56.406, 58.321, 60.656, 62.433, 65.078, 67.969, 71.239};
constexpr std::array<double, 59> kSemiMajor = {
    165.9, 166.7, 167.6, 168.4, 169.0, 169.0, 168.5, 167.8, 167.5, 167.8, 168.5, 169.0, 168.9, 168.2, 168.1,
    168.8, 169.3, 168.8, 167.4, 166.8, 167.5, 168.5, 168.2, 167.1, 166.5, 166.7, 166.0, 164.1, 162.9, 162.4,
    161.1, 159.3, 157.6, 155.4, 153.4, 151.7, 149.8, 148.1, 146.1, 143.4, 140.3, 137.0, 133.6, 130.8, 128.7,
    126.4, 124.0, 121.5, 118.4, 114.9, 111.6, 109.4, 107.2, 104.5, 100.8, 96.9,  93.0,  88.3,  82.7};
constexpr std::array<double, 59> kGlare = {
    28.6, 29.3, 29.8, 30.2, 30.3, 30.5, 30.9, 31.4, 32.0, 32.5, 32.9, 33.2, 33.6, 34.0, 34.4,
    34.7, 35.1, 35.7, 36.4, 36.8, 37.0, 37.3, 37.8, 38.4, 38.8, 39.1, 39.7, 40.4, 40.8, 41.1,
    41.6, 42.2, 42.6, 43.0, 43.4, 43.9, 44.3, 44.7, 45.2, 45.7, 46.2, 46.8, 47.3, 47.9, 48.5,
    49.0, 49.6, 50.2, 50.9, 51.6, 52.4, 53.1, 53.7, 54.2, 55.0, 55.9, 57.1, 58.9, 60.8};

double table(const std::array<double, 59>& v, double t) {
    const double x = std::clamp(t / kStep, 0.0, static_cast<double>(v.size() - 1));
    const size_t i = std::min(static_cast<size_t>(x), v.size() - 2);
    const double f = x - static_cast<double>(i);
    return v[i] + (v[i + 1] - v[i]) * f;
}

double smoothstep(double a, double b, double x) {
    const double u = std::clamp((x - a) / (b - a), 0.0, 1.0);
    return u * u * (3.0 - 2.0 * u);
}

// Pixels of the reference to lengths at the camera's distance: 1 frame height spans
// 1 / kFocal of them there.
double pixelsToLength(double px) {
    return px / (kFocal * kRefHeight);
}

// Where each shot's binary is centred (frame coordinates), fitted with its orbit.
constexpr double kOpeningX = -0.0140, kOpeningY = 0.0360;
constexpr double kClosingX = -0.0150, kClosingY = 0.0350;

glm::dvec3 centreAt(double x, double y) {
    return glm::dvec3(x / kFocal, y / kFocal, -1.0);
}

}  // namespace

glm::dvec3 orbitU() {
    return glm::dvec3(0.0, 1.0, 0.0);
}

glm::dvec3 orbitW() {
    return glm::dvec3(-kOrbitCos, 0.0, std::sqrt(1.0 - kOrbitCos * kOrbitCos));
}

Binary openingBinary(double t) {
    Binary b;
    b.centre = centreAt(kOpeningX, kOpeningY);
    const double end = kStep * static_cast<double>(kPhase.size() - 1);  // 5.8 s
    if (t <= end) {
        b.phase = table(kPhase, t);
        b.orbit = pixelsToLength(table(kSemiMajor, t));
        b.radius = pixelsToLength(table(kGlare, t));
    } else {
        // They touch: the two cores fuse into one spinning peanut, still whirling (30
        // radians a second, halving in 0.7 s), its lobes closing from 83 to 40 pixels apart
        // (fitted to the bright blob's second moments, 5.9 - 6.5 s).
        const double s = t - end;
        b.phase = kPhase.back() + 30.0 * 1.0 * (1.0 - std::exp(-s / 1.0));
        b.orbit = pixelsToLength(40.0 + (kSemiMajor.back() - 40.0) * std::exp(-s / 0.25));
        b.radius = pixelsToLength(kGlare.back() + 14.0 * smoothstep(0.0, 0.7, s));
    }
    b.merge = smoothstep(5.75, 6.1, t);
    return b;
}

Binary closingBinary(double t) {
    // The last shot, a wider view of the same tilted orbit, turning slowly: fitted to the
    // stars' positions from 15.6 to 16.2 s (semi-axis 153 pixels, glare 20).
    Binary b;
    b.centre = centreAt(kClosingX, kClosingY);
    b.phase = 1.56 + 1.35 * (t - 15.6);
    b.orbit = pixelsToLength(153.0);
    b.radius = pixelsToLength(20.1);
    return b;
}

glm::dvec3 starPosition(const Binary& b, int star) {
    const double sign = star == 0 ? 1.0 : -1.0;
    return b.centre + sign * b.orbit * (std::cos(b.phase) * orbitU() + std::sin(b.phase) * orbitW());
}

glm::dvec2 project(const glm::dvec3& p) {
    return glm::dvec2(p.x, p.y) * (kFocal / -p.z);
}

glm::dvec2 toReferencePixels(const glm::dvec2& frame) {
    return glm::dvec2(190.0 + frame.x * kRefHeight, 344.0 - frame.y * kRefHeight);
}

// --- Captions --------------------------------------------------------------------------------

const std::vector<Caption>& captions() {
    // Read off the reference frame by frame: when each appears and is replaced, its words and
    // their colours, how it comes in.
    static const std::vector<Caption> list = [] {
        using I = Ink;
        std::vector<Caption> c;
        Caption title;
        title.start = -1.0;
        title.end = 6.95;
        title.lines = {{{"This is what happens when", I::White}},
                       {{"two ", I::White}, {"neutron stars ", I::Yellow}, {"collide...", I::Red}}};
        title.leftAligned = true;
        title.refWidth = 335.0;
        title.entry = Entry::None;
        title.glitchOut = true;
        c.push_back(title);
        auto one = [&](double start, double end, double width, std::vector<Word> words, Entry entry, double settle) {
            Caption k;
            k.refWidth = width;
            k.start = start;
            k.end = end;
            k.lines = {std::move(words)};
            k.entry = entry;
            k.settle = settle;
            c.push_back(k);
        };
        one(7.15, 8.08, 311.0, {{"Emit ", I::Yellow}, {"gravitational waves", I::White}}, Entry::SlideFromLeft, 0.2);
        one(8.08, 9.02, 220.0, {{"Gamma-", I::White}, {"ray burst", I::Red}}, Entry::Glow, 0.25);
        one(9.02, 9.97, 183.0, {{"Merge ", I::Yellow}, {"into one", I::White}}, Entry::Glow, 0.4);
        one(10.08, 10.88, 289.0, {{"Create ", I::White}, {"heavy ", I::Yellow}, {"elements", I::White}}, Entry::Glow, 0.3);
        one(10.88, 11.95, 352.0, {{"Emit ", I::White}, {"electromagnetic ", I::Yellow}, {"radiation", I::Red}},
            Entry::SlideFromRight, 0.3);
        one(11.95, 12.88, 224.0, {{"Distort ", I::White}, {"spacetime", I::Red}}, Entry::Glow, 0.4);
        one(12.88, 13.82, 158.0, {{"Eject matter", I::White}}, Entry::Glow, 0.4);
        one(13.82, 14.72, 254.0, {{"Form new ", I::White}, {"black hole", I::Red}}, Entry::Glow, 0.4);
        Caption again = title;
        again.start = 14.88;
        again.end = kLength + 1.0;
        again.entry = Entry::SlideFromLeft;
        again.settle = 0.55;
        again.glitchOut = false;
        c.push_back(again);
        return c;
    }();
    return list;
}

const char* shotName(double t) {
    if (t < 5.85) return "two neutron stars spiral together";
    if (t < 6.75) return "they merge in a flash";
    if (t < 7.15) return "the merger's light";
    if (t < 8.08) return "gravitational waves ripple out; jets drill out";
    if (t < 9.02) return "the gamma-ray burst";
    if (t < 9.97) return "the debris merges into one cloud";
    if (t < 10.9) return "the cloud forges heavy elements";
    if (t < 11.95) return "the kilonova glows";
    if (t < 12.88) return "spacetime rings";
    if (t < 13.82) return "matter is ejected";
    if (t < 14.75) return "a new black hole";
    return "two neutron stars, before it all";
}

}  // namespace film
