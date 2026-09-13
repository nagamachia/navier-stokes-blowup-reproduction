#include "appendix_a7_log_discrepancy.hpp"

#include <cmath>
#include <iostream>

int main() {
    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;

    const auto q = nsblowup::appendix_a7_leading_heat_log_discrepancy(p, 0.4);
    if (!(q.Cp.sign == -1 && q.S.sign == 1 && q.I.sign == -1)) {
        std::cerr << "unexpected leading heat discrepancy signs\n";
        return 1;
    }
    if (!(q.patch_target_Cp.sign == 1 && q.patch_target_S.sign == -1 &&
          q.patch_target_I.sign == 1)) {
        std::cerr << "unexpected compensation target signs\n";
        return 2;
    }
    if (!(std::isfinite(q.patch_target_Cp.log_abs) &&
          std::isfinite(q.patch_target_S.log_abs) &&
          std::isfinite(q.patch_target_I.log_abs))) {
        std::cerr << "non-finite A.7 patch target logs\n";
        return 3;
    }
    if (!(q.patch_target_Cp.log_abs < -700.0 &&
          q.patch_target_S.log_abs < -700.0 &&
          q.patch_target_I.log_abs < -700.0)) {
        std::cerr << "ordered A.7 compensation target should be below double range\n";
        return 4;
    }

    const auto err = nsblowup::appendix_a7_leading_heat_error_audit(p, 0.4);
    if (!(std::isfinite(err.log_relative_bound_CpS) &&
          std::isfinite(err.log_relative_bound_I))) {
        std::cerr << "non-finite A.7 leading heat error bound\n";
        return 5;
    }
    if (!(err.log_relative_bound_CpS < -1.0e6 &&
          err.log_relative_bound_I < -1.0e6)) {
        std::cerr << "A.7 leading heat expansion is not sufficiently accurate in ordered regime\n";
        return 6;
    }

    return 0;
}
