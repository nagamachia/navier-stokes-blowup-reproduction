#pragma once

#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct CartesianLeadingField {
    double u1;
    double u2;
    double u3;
    double pressure;
};

// Equations (4.4)-(4.5), written in the axis-regular variables
// E=sqrt(2X)F and V0=X v0, with X=(x1^2+x2^2)/(2q).
template <class FFunc, class UFunc, class V0RegularFunc, class PiFunc>
CartesianLeadingField cartesian_leading_field(
    FFunc&& F,
    UFunc&& U,
    V0RegularFunc&& v0,
    PiFunc&& Pi,
    double x1,
    double x2,
    double q,
    double eta,
    double h) {
    if (!(q > 0.0)) {
        throw std::invalid_argument("q must be positive");
    }
    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }

    const double A = 0.5 + h;
    const double X = (x1 * x1 + x2 * x2) / (2.0 * q);
    const double radial_coefficient = v0(X, eta) / (2.0 * q);
    const double swirl_coefficient = std::pow(q, -A - 0.5) * F(X, eta);

    return {
        radial_coefficient * x1 - swirl_coefficient * x2,
        radial_coefficient * x2 + swirl_coefficient * x1,
        std::pow(q, -A) * U(X, eta),
        std::pow(q, -2.0 * A) * Pi(X, eta),
    };
}

}  // namespace nsblowup
