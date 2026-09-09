#pragma once

#include <stdexcept>

namespace nsblowup {

struct ProfileJet {
    double value;
    double d_eta;
    double d_x_dilation;  // X * partial_X f
};

inline double T_b(double b, double h, double eta, const ProfileJet& jet) {
    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }
    const double D = 0.5 - h;
    const double L = 1.0 - 2.0 * h * eta * eta;
    return (-b * jet.value
            + D * eta * jet.d_eta
            + jet.d_x_dilation) / L;
}

inline double Z_b(double b, double h, double eta, const ProfileJet& jet) {
    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }
    const double d = 1.0 - eta * eta;
    const double L = 1.0 - 2.0 * h * eta * eta;
    return (2.0 * b * eta * jet.value
            + d * jet.d_eta
            - 2.0 * eta * jet.d_x_dilation) / L;
}

}  // namespace nsblowup
