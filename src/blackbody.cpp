#include "blackbody.h"

#include <algorithm>
#include <cmath>

namespace blackbody {

namespace {

// Piecewise Gaussian used by the CIE fit below.
double lobe(double x, double mu, double s1, double s2) {
    double t = (x - mu) / (x < mu ? s1 : s2);
    return std::exp(-0.5 * t * t);
}

// Multi-lobe fit of the CIE 1931 2-degree observer
// (Wyman, Sloan & Shirley 2013), lambda in nanometres.
glm::dvec3 colorMatch(double lambda) {
    return {1.056 * lobe(lambda, 599.8, 37.9, 31.0) + 0.362 * lobe(lambda, 442.0, 16.0, 26.7) -
                0.065 * lobe(lambda, 501.1, 20.4, 26.2),
            0.821 * lobe(lambda, 568.8, 46.9, 40.5) + 0.286 * lobe(lambda, 530.9, 16.3, 31.1),
            1.217 * lobe(lambda, 437.0, 11.8, 36.0) + 0.681 * lobe(lambda, 459.0, 26.0, 13.8)};
}

glm::dvec3 planckXYZ(double temperature) {
    constexpr double h = 6.62607015e-34, c = 2.99792458e8, k = 1.380649e-23;
    glm::dvec3 xyz(0.0);
    for (double nm = 360.0; nm <= 830.0; nm += 2.0) {
        double lambda = nm * 1e-9;
        double radiance = 2.0 * h * c * c / std::pow(lambda, 5.0) / std::expm1(h * c / (lambda * k * temperature));
        xyz += radiance * colorMatch(nm);
    }
    return xyz;
}

glm::dvec3 xyzToLinearSrgb(const glm::dvec3& v) {
    glm::dvec3 rgb(3.2406 * v.x - 1.5372 * v.y - 0.4986 * v.z,
                   -0.9689 * v.x + 1.8758 * v.y + 0.0415 * v.z,
                   0.0557 * v.x - 0.2040 * v.y + 1.0570 * v.z);
    return glm::max(rgb, glm::dvec3(0.0));  // out-of-gamut reds have small negative components
}

}  // namespace

std::vector<float> makeTable(int n) {
    const double reference = planckXYZ(6500.0).y;
    std::vector<float> table;
    table.reserve(static_cast<size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
        double logT = kLogTMin + (kLogTMax - kLogTMin) * i / (n - 1);
        glm::dvec3 rgb = xyzToLinearSrgb(planckXYZ(std::pow(10.0, logT))) / reference;
        table.push_back(static_cast<float>(rgb.r));
        table.push_back(static_cast<float>(rgb.g));
        table.push_back(static_cast<float>(rgb.b));
    }
    return table;
}

double luminance(double temperature) {
    return planckXYZ(temperature).y / planckXYZ(6500.0).y;
}

glm::vec3 hue(double temperature) {
    glm::dvec3 rgb = xyzToLinearSrgb(planckXYZ(temperature));
    double peak = std::max({rgb.r, rgb.g, rgb.b, 1e-300});
    return glm::vec3(rgb / peak);
}

}  // namespace blackbody
