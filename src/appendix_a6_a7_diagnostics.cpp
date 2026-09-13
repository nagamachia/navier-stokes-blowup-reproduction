#include "appendix_a6_heat.hpp"
#include "appendix_a7_heat_compensation.hpp"
#include "appendix_a7_log_compensation.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    using namespace nsblowup;
    const std::string out = argc > 1 ? argv[1] : "appendix_a6_a7_heat.csv";
    std::ofstream os(out);
    if (!os) {
        std::cerr << "failed to open " << out << "\n";
        return 1;
    }

    os << "eta,heat_ratio_x100,heat_ode_residual,jacobian_det_scaled,"
          "recovery_error_inf,c0,c1,c2,ordered_target_log_cp,"
          "ordered_target_log_s,ordered_target_log_i,target_log_spread,"
          "quadratic_to_dominant_log_ratio,quadratic_log_magnitude,"
          "leading_error_log_bound_cps,leading_error_log_bound_i,"
          "all_components_resolvable_in_double\n";
    os << std::setprecision(15);

    constexpr double h = 0.08;
    constexpr double lambda = 0.03;
    const Vector known{1.2e-4, -8.0e-5, 5.0e-5};

    OuterProfileParameters ordered;
    ordered.lambda = 1e-5;
    ordered.log_h = -80.0;
    ordered.log_Pstar = 70.0;

    for (int k = -9; k <= 9; ++k) {
        const double eta = static_cast<double>(k) / 10.0;
        const double Z = 2.0 * (1.0 - eta * eta) / 100.0;
        const double heat_ratio = appendix_a6_heat_ratio(100.0, eta, h);
        const double ode = appendix_a6_heat_ode_residual(Z, h);

        AppendixA7HeatCompensation patch(eta, lambda);
        Matrix B = patch.jacobian_at_zero();
        for (std::size_t i = 0; i < B.size(); ++i) {
            double scale = 0.0;
            for (double v : B[i]) scale = std::max(scale, std::abs(v));
            for (double& v : B[i]) v /= scale;
        }
        const double det = determinant(B);
        const Vector target = patch.change(known).vector();
        const Vector solved = patch.solve_for_change(target);
        const Vector recovered = patch.change(solved).vector();
        double err = 0.0;
        for (std::size_t i = 0; i < target.size(); ++i)
            err = std::max(err, std::abs(recovered[i] - target[i]));

        const auto log_target = appendix_a7_leading_heat_log_discrepancy(ordered, eta);
        const auto hierarchy = appendix_a7_audit_ordered_heat_compensation(
            ordered, eta, 20.0, 0.05, 1024);
        const auto asymptotic = appendix_a7_leading_heat_error_audit(ordered, eta);

        os << eta << ',' << heat_ratio << ',' << ode << ',' << det << ','
           << err << ',' << solved[0] << ',' << solved[1] << ',' << solved[2] << ','
           << log_target.patch_target_Cp.log_abs << ','
           << log_target.patch_target_S.log_abs << ','
           << log_target.patch_target_I.log_abs << ','
           << hierarchy.target_log_spread << ','
           << hierarchy.quadratic_to_dominant_log_ratio << ','
           << hierarchy.quadratic_log_magnitude << ','
           << asymptotic.log_relative_bound_CpS << ','
           << asymptotic.log_relative_bound_I << ','
           << (hierarchy.all_components_resolvable_in_double ? 1 : 0) << '\n';
    }

    std::cout << "wrote Appendix A.6/A.7 diagnostics to " << out << "\n";
    return 0;
}
