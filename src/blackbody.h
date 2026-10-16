#pragma once

#include <glm/glm.hpp>

#include <vector>

// Colors of thermal (Planck) radiation, integrated against the CIE 1931 color
// matching functions and converted to linear sRGB.
namespace blackbody {

constexpr double kLogTMin = 2.5;  // log10 of the lowest tabulated temperature (316 K)
constexpr double kLogTMax = 5.3;  // log10 of the highest (200,000 K)

// n RGB triples for temperatures spaced evenly in log10(T) from kLogTMin to kLogTMax.
// Values are absolute visible-band radiance, scaled so a 6500 K blackbody has
// luminance 1. Hotter emitters are therefore brighter as well as bluer.
std::vector<float> makeTable(int n);

// Visible luminance of a blackbody at temperature T, relative to one at 6500 K.
double luminance(double temperature);

// Hue of a blackbody at temperature T, normalized so the largest channel is 1.
glm::vec3 hue(double temperature);

}  // namespace blackbody
