#include "tidal.h"

#include <algorithm>
#include <cmath>

namespace tidal {

Mat3 identity() {
    Mat3 m{};
    for (int i = 0; i < 3; ++i) m[i][i] = 1.0;
    return m;
}

Mat3 multiply(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) out[i][j] += a[i][k] * b[k][j];
        }
    }
    return out;
}

Mat3 transpose(const Mat3& m) {
    Mat3 out{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) out[i][j] = m[j][i];
    }
    return out;
}

double determinant(const Mat3& m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

Mat3 inverse(const Mat3& m) {
    const double inv = 1.0 / determinant(m);
    Mat3 out;
    out[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv;
    out[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv;
    out[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv;
    out[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv;
    out[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv;
    out[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv;
    out[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv;
    out[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv;
    out[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv;
    return out;
}

double carlsonRD(double x, double y, double z) {
    // Duplication algorithm (Carlson 1995).
    double sum = 0.0, factor = 1.0;
    for (int i = 0; i < 60; ++i) {
        const double mu = (x + y + 3.0 * z) / 5.0;
        const double dx = (mu - x) / mu, dy = (mu - y) / mu, dz = (mu - z) / mu;
        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) < 1e-4) {
            const double ea = dx * dy, eb = dz * dz, ec = ea - eb, ed = ea - 6.0 * eb, ee = ed + ec + ec;
            const double series = 1.0 + ed * (-3.0 / 14.0 + 9.0 / 88.0 * ed - 4.5 / 26.0 * dz * ee) +
                                  dz * (ee / 6.0 + dz * (-9.0 / 22.0 * ec + dz * 3.0 / 26.0 * ea));
            return 3.0 * sum + factor * series / (mu * std::sqrt(mu));
        }
        const double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
        const double lambda = sx * (sy + sz) + sy * sz;
        sum += factor / (sz * (z + lambda));
        factor *= 0.25;
        x = 0.25 * (x + lambda);
        y = 0.25 * (y + lambda);
        z = 0.25 * (z + lambda);
    }
    const double mu = (x + y + 3.0 * z) / 5.0;
    return 3.0 * sum + factor / (mu * std::sqrt(mu));
}

void symmetricEigen(const Mat3& m, Vec3& values, Mat3& vectors) {
    // Cyclic Jacobi rotations.
    Mat3 a = m;
    vectors = identity();
    for (int sweep = 0; sweep < 50; ++sweep) {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        const double scale = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2];
        if (off <= 1e-30 * std::max(scale, 1e-300)) break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(a[p][q]) < 1e-300) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    const double vkp = vectors[k][p], vkq = vectors[k][q];
                    vectors[k][p] = c * vkp - s * vkq;
                    vectors[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    // Sort ascending.
    int order[3] = {0, 1, 2};
    std::sort(order, order + 3, [&](int i, int j) { return a[i][i] < a[j][j]; });
    Mat3 sorted;
    for (int j = 0; j < 3; ++j) {
        values[j] = a[order[j]][order[j]];
        for (int k = 0; k < 3; ++k) sorted[k][j] = vectors[k][order[j]];
    }
    vectors = sorted;
}

Vec3 indexSymbols(const Vec3& a) {
    // A_i = a1 a2 a3 integral_0^inf du / ((a_i^2 + u) Delta) = (2/3) a1 a2 a3 R_D(a_j^2, a_k^2, a_i^2).
    const double volume = a[0] * a[1] * a[2];
    const double s[3] = {a[0] * a[0], a[1] * a[1], a[2] * a[2]};
    return {2.0 / 3.0 * volume * carlsonRD(s[1], s[2], s[0]), 2.0 / 3.0 * volume * carlsonRD(s[0], s[2], s[1]),
            2.0 / 3.0 * volume * carlsonRD(s[0], s[1], s[2])};
}

Mat3 tidalTensor(double a, const kerr::Constants& c, const kerr::State& s) {
    const auto x = kerr::toCartesian(a, s.r, s.theta, s.phi);
    const kerr::Derivatives d = kerr::derivatives(a, c, s);
    const auto v = kerr::cartesianVelocity(a, s.r, s.theta, s.phi, d.dr, d.dtheta, d.dphi);
    const double xl = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    const Vec3 n = {x[0] / xl, x[1] / xl, x[2] / xl};
    Vec3 l = {n[1] * v[2] - n[2] * v[1], n[2] * v[0] - n[0] * v[2], n[0] * v[1] - n[1] * v[0]};
    double ll = std::sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    if (ll < 1e-9) {
        // Radial motion: any normal will do.
        const Vec3 z = std::abs(n[2]) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
        l = {n[1] * z[2] - n[2] * z[1], n[2] * z[0] - n[0] * z[2], n[0] * z[1] - n[1] * z[0]};
        ll = std::sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    }
    for (double& e : l) e /= ll;
    const Vec3 m = {l[1] * n[2] - l[2] * n[1], l[2] * n[0] - l[0] * n[2], l[0] * n[1] - l[1] * n[0]};

    const double K = c.Q + (c.L - a * c.E) * (c.L - a * c.E);
    const double k = 3.0 * K / (s.r * s.r);
    const double unit = 1.0 / (s.r * s.r * s.r);
    Mat3 C{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = unit * ((2.0 + k) * n[i] * n[j] - (1.0 + k) * m[i] * m[j] - l[i] * l[j]);
        }
    }
    return C;
}

Mat3 acceleration(const Mat3& F, const Mat3& C, double w2, double gamma) {
    const double det = determinant(F);
    // Self-gravity: -(3/2) w^2 / det F * A F, with A from the ellipsoid's axes.
    const Mat3 B = multiply(F, transpose(F));
    Vec3 lambda;
    Mat3 R;
    symmetricEigen(B, lambda, R);
    const Vec3 axes = {std::sqrt(std::max(lambda[0], 1e-12)), std::sqrt(std::max(lambda[1], 1e-12)),
                       std::sqrt(std::max(lambda[2], 1e-12))};
    const Vec3 A = indexSymbols(axes);
    Mat3 Am{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) Am[i][j] += R[i][k] * A[k] * R[j][k];
        }
    }
    const Mat3 CF = multiply(C, F), AF = multiply(Am, F), invT = transpose(inverse(F));
    const double gravity = 1.5 * w2 / det, pressure = w2 * std::pow(det, 1.0 - gamma);
    Mat3 out;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) out[i][j] = CF[i][j] - gravity * AF[i][j] + pressure * invT[i][j];
    }
    return out;
}

void advance(Star& star, const Mat3& C, double w2, double dtau) {
    // Steps well inside the star's own oscillation period, 2 pi / w.
    const int n = std::max(1, static_cast<int>(std::ceil(std::abs(dtau) * std::sqrt(w2) / 0.02)));
    const double h = dtau / n;
    Mat3 acc = acceleration(star.F, C, w2);
    for (int step = 0; step < n; ++step) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                star.Fdot[i][j] += 0.5 * h * acc[i][j];
                star.F[i][j] += h * star.Fdot[i][j];
            }
        }
        acc = acceleration(star.F, C, w2);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) star.Fdot[i][j] += 0.5 * h * acc[i][j];
        }
    }
}

double stretch(const Mat3& F) {
    Vec3 lambda;
    Mat3 R;
    symmetricEigen(multiply(F, transpose(F)), lambda, R);
    return std::sqrt(std::max(lambda[2], 0.0));
}

Mat3 inverseShape(const Mat3& F) {
    Vec3 lambda;
    Mat3 R;
    symmetricEigen(multiply(F, transpose(F)), lambda, R);
    Mat3 out{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) out[i][j] += R[i][k] * R[j][k] / std::sqrt(std::max(lambda[k], 1e-12));
        }
    }
    return out;
}

}  // namespace tidal
