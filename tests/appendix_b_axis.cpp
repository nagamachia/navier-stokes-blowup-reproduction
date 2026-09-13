#include "appendix_b_axis_asymptotic.hpp"

#include <cmath>
#include <iostream>

int main() {
    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;

    const auto sep = nsblowup::appendix_b_separation_audit(p, 0.05, 20.0, 0.05, 0.04, 161);
    if (!sep.unique_H_zero || !(sep.eta0 < 0.0 && sep.eta0 > -0.1)) {
        std::cerr << "B.1 H* zero is not in the expected interval\n";
        return 1;
    }
    if (!sep.positive_Z_at_H_zero || !(sep.normalized_Zstar_eta0 > 0.0)) {
        std::cerr << "B.1 Z*(eta0) positivity failed\n";
        return 2;
    }
    if (!sep.chi_separation_certified || !(sep.min_chi_on_small_Z > 0.99)) {
        std::cerr << "B.2 chi separation failed\n";
        return 3;
    }

    for (double eta : {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto a = nsblowup::appendix_b_axis_data(p, eta);
        if (!std::isfinite(a.normalized_Zstar)) {
            std::cerr << "B.1 normalized Z* is not finite\n";
            return 4;
        }
        if (!(-a.Wstar > 2.8)) {
            std::cerr << "B.1 -W*>2.8 bound failed\n";
            return 5;
        }
        const double log_phi = nsblowup::appendix_b_log_phi_star(
            p, eta, 25.0, sep.sigma_star, 0.05, 400);
        if (!std::isfinite(log_phi)) {
            std::cerr << "B.3 log phi* is not finite\n";
            return 6;
        }
    }

    double f0min = 1.0;
    for (int i = 0; i <= 410; ++i) {
        const double z = 4.1 * static_cast<double>(i) / 410.0;
        const double f0 = nsblowup::appendix_b_f0(z);
        f0min = std::min(f0min, f0);
        if (!(f0 > 0.265 && f0 <= 1.0 + 1e-14)) {
            std::cerr << "B.11 comparison profile bound failed at z=" << z << '\n';
            return 7;
        }
    }
    if (!(f0min > 0.265)) return 8;

    const auto endpoint = nsblowup::appendix_b_endpoint_asymptotic_audit(
        p, sep, 25.0, 0.05, 121);
    if (!(endpoint.min_phi0 > 0.265)) {
        std::cerr << "B.3 comparison profile lost positivity\n";
        return 9;
    }
    if (!endpoint.azimuthal_branch_certified ||
        !(endpoint.min_azimuthal_p1_on_chi_branch > 2.36)) {
        std::cerr << "B.19 azimuthal endpoint branch failed\n";
        return 10;
    }
    if (!endpoint.axial_branch_separated ||
        !(endpoint.min_abs_normalized_ns_on_axial_branch > sep.normalized_delta_star)) {
        std::cerr << "B.19 axial endpoint branch failed separation\n";
        return 11;
    }
    if (!std::isfinite(endpoint.max_log_required_C_normalized)) {
        std::cerr << "B.19 normalized amplitude threshold is not finite\n";
        return 12;
    }

    return 0;
}
