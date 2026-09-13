#include "appendix_a6_heat.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double atol, double rtol) {
    return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}

}

int main() {
    using namespace nsblowup;

    const double h = 0.08;
    const double H0 = appendix_a6_heat_factor(0.0, h);
    if (!close(H0, 1.0, 3e-5, 3e-5)) {
        std::cerr << "A.32 normalization failed: H(0)=" << H0 << "\n";
        return 1;
    }

    for (int m = 1; m <= 2; ++m) {
        const double numeric = appendix_a6_heat_factor_derivative(0.0, h, m);
        const double exact = appendix_a6_heat_factor_at_zero_derivative(h, m);
        if (!close(numeric, exact, 8e-5, 8e-5)) {
            std::cerr << "A.35 derivative mismatch m=" << m
                      << " numeric=" << numeric << " exact=" << exact << "\n";
            return 1;
        }
    }

    for (double Z : {0.02, 0.2, 0.8, 2.0}) {
        const double H = appendix_a6_heat_factor(Z, h);
        const double Hp = appendix_a6_heat_factor_derivative(Z, h, 1);
        const double residual = appendix_a6_heat_ode_residual(Z, h);
        if (!(H > 0.0) || !(Hp < 0.0)) {
            std::cerr << "A.6 positivity/monotonicity failed at Z=" << Z << "\n";
            return 1;
        }
        if (std::abs(residual) > 2e-4) {
            std::cerr << "A.37 ODE residual too large at Z=" << Z
                      << ": " << residual << "\n";
            return 1;
        }
        const double log_slope = -Z * Hp / H;
        if (!(log_slope >= 0.0 && log_slope < h + 2e-4)) {
            std::cerr << "A.6 logarithmic slope bound failed at Z=" << Z << "\n";
            return 1;
        }
    }

    const double eta = 0.4;
    const double d = 1.0 - eta * eta;
    for (double X : {80.0, 160.0, 320.0}) {
        const double ratio = appendix_a6_heat_ratio(X, eta, h);
        const double first = 1.0 - 2.0 * h * (1.0 + h) * d / X;
        if (!close(ratio, first, 8e-6, 8e-6)) {
            std::cerr << "A.38 expansion mismatch at X=" << X
                      << ": ratio=" << ratio << " first=" << first << "\n";
            return 1;
        }
    }

    const double X = 50.0;
    const double heat = appendix_a6_heat_ratio(X, eta, h);
    if (!close(appendix_a6_terminal_replacement_ratio(X, eta, h, 0.1),
               1.0, 1e-12, 1e-12)) {
        std::cerr << "A.39 changed the profile before cutoff\n";
        return 1;
    }
    if (!close(appendix_a6_terminal_replacement_ratio(X, eta, h, 0.6),
               heat, 1e-12, 1e-12)) {
        std::cerr << "A.39 failed to reach full heat profile\n";
        return 1;
    }

    if (!close(appendix_a6_heat_ratio(20.0, 1.0, h), H0, 1e-10, 1e-10)) {
        std::cerr << "eta endpoint regularity failed\n";
        return 1;
    }

    std::cout << "Appendix A.6 heat regression passed\n";
    return 0;
}
