#include "appendix_a7_log_compensation.hpp"

#include <cmath>
#include <iostream>

int main() {
    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;

    const auto q = nsblowup::appendix_a7_leading_heat_log_discrepancy(p, 0.4);
    const auto sol = nsblowup::appendix_a7_solve_ordered_heat_compensation(p, 0.4, 20.0, 0.05, 1024);

    if (!(std::isfinite(sol.common_log_scale) && sol.common_log_scale < -700.0)) {
        std::cerr << "A.7 compensation common scale should be below double range\n";
        return 1;
    }
    if (nsblowup::norm_inf(sol.scaled_linear_residual) > 1e-10) {
        std::cerr << "A.7 scaled linear compensation residual too large\n";
        return 2;
    }
    for (const auto& c : sol.coefficients) {
        if (c.sign != 0 && !std::isfinite(c.log_abs)) {
            std::cerr << "non-finite A.7 signed-log coefficient\n";
            return 3;
        }
    }
    if (!(sol.quadratic_to_linear_log_bound < -650.0)) {
        std::cerr << "A.7 quadratic correction is not negligible on ordered scale\n";
        return 4;
    }

    // Cp is the axis pressure datum increment, so solving the same three-vector
    // also restores the pressure datum at leading order in the common scale.
    const double cp_scaled = static_cast<double>(q.patch_target_Cp.sign) *
        std::exp(q.patch_target_Cp.log_abs - sol.common_log_scale);
    if (std::abs(sol.scaled_target[0] - cp_scaled) > 1e-14) {
        std::cerr << "A.7 pressure target was not passed through consistently\n";
        return 5;
    }

    return 0;
}
