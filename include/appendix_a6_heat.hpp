#pragma once

#include "outer_profile_stages.hpp"

#include <cmath>
#include <stdexcept>

namespace nsblowup {

// Appendix A.6, equations (A.32)--(A.39).
inline double appendix_a6_rising_factorial(double b, int m) {
    double out = 1.0;
    for (int k = 0; k < m; ++k) out *= (b + static_cast<double>(k));
    return out;
}

inline double appendix_a6_heat_factor_derivative(double Z, double h, int m,
                                                  int panels = 8192,
                                                  double vmax = 48.0) {
    if (!(Z >= 0.0)) throw std::invalid_argument("A.32 requires Z >= 0");
    if (!(h > 0.0 && h < 0.5)) throw std::invalid_argument("A.6 requires 0 < h < 1/2");
    if (m < 0) throw std::invalid_argument("derivative order must be nonnegative");
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    const double prefactor = ((m % 2) ? -1.0 : 1.0) *
        appendix_a6_rising_factorial(h, m) / std::tgamma(1.0 + h);
    const double dv = vmax / static_cast<double>(panels);
    auto integrand = [&](double v) {
        if (v == 0.0) return 0.0;
        return std::exp(-v) * std::pow(v, h + static_cast<double>(m)) *
               std::pow(1.0 + Z * v, -h - static_cast<double>(m));
    };
    double sum = 0.0;
    for (int i = 0; i <= panels; ++i) {
        const double v = dv * static_cast<double>(i);
        const double w = (i == 0 || i == panels) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum += w * integrand(v);
    }
    return prefactor * sum * dv / 3.0;
}

inline double appendix_a6_heat_factor(double Z, double h) {
    return appendix_a6_heat_factor_derivative(Z, h, 0);
}

inline double appendix_a6_heat_factor_at_zero_derivative(double h, int m) {
    if (!(h > 0.0 && h < 0.5)) throw std::invalid_argument("A.6 requires 0 < h < 1/2");
    if (m < 0) throw std::invalid_argument("derivative order must be nonnegative");
    return ((m % 2) ? -1.0 : 1.0) * appendix_a6_rising_factorial(h, m) *
           appendix_a6_rising_factorial(1.0 + h, m);
}

inline double appendix_a6_heat_ratio(double X, double eta, double h) {
    if (!(X > 0.0)) throw std::invalid_argument("X must be positive");
    if (std::abs(eta) > 1.0) throw std::invalid_argument("|eta| must not exceed 1");
    return appendix_a6_heat_factor(2.0 * (1.0 - eta * eta) / X, h);
}

inline double appendix_a6_terminal_replacement_ratio(double X, double eta,
                                                       double h, double y) {
    const double chi = appendix_a_smooth_step((y - 0.2) / 0.3);
    return 1.0 + chi * (appendix_a6_heat_ratio(X, eta, h) - 1.0);
}

inline double appendix_a6_heat_ode_residual(double Z, double h) {
    const double aK = 1.0 + h;
    const double H = appendix_a6_heat_factor_derivative(Z, h, 0);
    const double Hp = appendix_a6_heat_factor_derivative(Z, h, 1);
    const double Hpp = appendix_a6_heat_factor_derivative(Z, h, 2);
    return Z * Z * Hpp + (1.0 + 2.0 * aK * Z) * Hp +
           aK * (aK - 1.0) * H;
}

}  // namespace nsblowup
