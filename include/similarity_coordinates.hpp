#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct SimilarityCoordinates {
    double q;
    double eta;
    double X;
};

inline double solve_q(double tau, double z, double h) {
    if (!(tau > 0.0)) {
        throw std::invalid_argument("tau must be positive");
    }
    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }

    const double z2 = z * z;
    if (z2 == 0.0) {
        return tau;
    }

    // From tau = q(1-eta^2), z = q^(1/2-h) eta:
    //     q - z^2 q^(2h) = tau.
    // On the physical branch |eta| < 1 the derivative is
    // 1 - 2 h eta^2 >= 1 - 2h > 0, so the root is unique.
    auto f = [=](double q) {
        return q - z2 * std::pow(q, 2.0 * h) - tau;
    };

    double lo = tau;
    double hi = std::max({1.0, 2.0 * tau, tau + z2 + 1.0});
    while (f(hi) <= 0.0) {
        hi *= 2.0;
        if (!std::isfinite(hi)) {
            throw std::runtime_error("failed to bracket similarity-coordinate root");
        }
    }

    double q = 0.5 * (lo + hi);
    for (int iter = 0; iter < 120; ++iter) {
        q = 0.5 * (lo + hi);
        if (f(q) > 0.0) {
            hi = q;
        } else {
            lo = q;
        }
        if ((hi - lo) <= 2e-14 * q) {
            break;
        }
    }

    return 0.5 * (lo + hi);
}

inline SimilarityCoordinates coordinates(
    double r,
    double z,
    double tau,
    double h) {
    const double q = solve_q(tau, z, h);
    const double D = 0.5 - h;
    const double eta = z / std::pow(q, D);
    const double X = r * r / (2.0 * q);
    return {q, eta, X};
}

}  // namespace nsblowup
