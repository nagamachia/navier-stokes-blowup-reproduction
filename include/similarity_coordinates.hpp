#pragma once

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

    auto f = [=](double q) {
        return q - z2 / std::pow(q, 2.0 * h) - tau;
    };

    double lo = tau;
    double hi = tau + z2 / std::pow(tau, 2.0 * h) + 1.0;
    double q = 0.5 * (lo + hi);

    for (int iter = 0; iter < 100; ++iter) {
        const double value = f(q);
        if (value > 0.0) {
            hi = q;
        } else {
            lo = q;
        }

        const double derivative =
            1.0 + 2.0 * h * z2 / std::pow(q, 2.0 * h + 1.0);
        const double newton = q - value / derivative;
        const double midpoint = 0.5 * (lo + hi);
        q = (newton > lo && newton < hi) ? newton : midpoint;

        if ((hi - lo) <= 2e-14 * q) {
            break;
        }
    }

    return q;
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
