#pragma once

// The neutron star merger as the reference film shows it: the merger's "Film" view.
//
// The reference is a 17.6-second portrait (9:16) short, `refrence photos/neytron star
// collision recording/Screen Recording 2026-10-19 225954.mp4`, cut from NASA Goddard's
// animations of GW170817: two white-hot neutron stars whirl together trailing ribbons of
// light, merge in a flash, a jet bursts out between rings of gravitational waves, the
// debris billows into a cloud, then a glowing shell turns from orange to blue. The Film view
// replays those shots with their timing, framing and colours, measured frame by frame
// (shaders/merger_film.frag draws them); unlike the cinematic and simulation views it does
// not follow the hydrodynamics. The orbit of the first and last shots is fitted to the
// stars' tracked positions, and tests/test_kerr.cpp checks it against them.
//
// Frame coordinates: the 9:16 frame is 1 high, x to the right and y up from its centre;
// the reference's 380 x 688 pixels map to x = (px - 190) / 688, y = (344 - py) / 688.
// The camera sits at the origin looking down -z with a 45-degree vertical field of view.

#include <glm/glm.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace film {

constexpr double kLength = 17.6;        // s: the reference's length; the film loops
constexpr double kAspect = 380.0 / 688.0;  // its frame, width over height
constexpr double kFocal = 1.2071067811865475;  // 1 / (2 tan 22.5 deg), in frame heights
constexpr double kRefHeight = 688.0;     // the reference's frame height, pixels

// The tilt of the stars' orbit: seen 65 degrees from face-on, its circle is drawn as an
// ellipse 0.42 times as wide as it is tall (the ratio of the stars' narrowest and widest
// projected separations), the tall axis upright and the near side on the left.
constexpr double kOrbitCos = 0.42;

// The binary at a moment of the first or last shot. Star 1 is at
// centre + orbit (cos phase U + sin phase W), star 2 opposite, with U up and W = (-0.42, 0,
// 0.91) in camera space; lengths in units of the camera's distance to the centre.
struct Binary {
    glm::dvec3 centre{0.0, 0.0, -1.0};
    double phase = 0.0;
    double orbit = 0.2;     // the orbit's radius
    double radius = 0.035;  // how far a star's white glare reaches (its disk is a little smaller)
    double merge = 0.0;     // 0 apart, 1 fused into one
};
Binary openingBinary(double t);  // the first shot, 0 to 6.6 s
Binary closingBinary(double t);  // the last, 14.9 to 17.6 s
glm::dvec3 orbitU();
glm::dvec3 orbitW();
glm::dvec3 starPosition(const Binary& b, int star);  // camera space; star 0 or 1
// Where a camera-space point is drawn, in frame coordinates, and in the reference's pixels.
glm::dvec2 project(const glm::dvec3& p);
glm::dvec2 toReferencePixels(const glm::dvec2& frame);

// --- Captions ------------------------------------------------------------------------------

enum class Ink { White, Yellow, Red };
struct Word {
    std::string text;
    Ink ink = Ink::White;
};
// How a caption comes in: blurred and glowing with its letters spread, or sliding in from
// one side, settling over `settle` seconds.
enum class Entry { Glow, SlideFromLeft, SlideFromRight, None };
struct Caption {
    double start = 0.0, end = 0.0;  // s
    std::vector<std::vector<Word>> lines;
    bool leftAligned = false;       // the two-line title is set flush left
    double refWidth = 0.0;          // its widest line's width in the reference, pixels
    Entry entry = Entry::Glow;
    double settle = 0.35;
    bool glitchOut = false;         // leaves by smearing, as the first title does at the merger
};
const std::vector<Caption>& captions();

// Which shot is showing (for the panel and the logs).
const char* shotName(double t);

}  // namespace film
