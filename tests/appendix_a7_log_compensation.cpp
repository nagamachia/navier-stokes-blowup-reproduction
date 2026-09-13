#include "appendix_a7_log_compensation.hpp"

#include <cmath>
#include <iostream>

int main() {
    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;

    const auto audit = nsblowup::appendix_a7_audit_ordered_heat_compensation(
        p, 0.4, 20.0, 0.05, 1024);

    if (!(std::isfinite(audit.largest_target_log) &&
          audit.largest_target_log < -700.0)) {
        std::cerr << "A.7 dominant compensation scale should be below double range\n";
        return 1;
    }
    if (!(audit.target_log_spread > 1.0e6)) {
        std::cerr << "A.7 ordered target hierarchy is unexpectedly narrow\n";
        return 2;
    }
    if (audit.all_components_resolvable_in_double) {
        std::cerr << "A.7 ordered targets must not be labeled double-resolvable\n";
        return 3;
    }
    if (nsblowup::norm_inf(audit.dominant_scaled_linear_residual) > 1e-10) {
        std::cerr << "A.7 dominant scaled linear residual too large\n";
        return 4;
    }
    if (!(audit.quadratic_to_dominant_log_ratio < -650.0)) {
        std::cerr << "A.7 quadratic term is not negligible relative to dominant target\n";
        return 5;
    }
    if (!(audit.quadratic_log_magnitude > audit.smallest_target_log)) {
        std::cerr << "A.7 hierarchy audit expected quadratic term to exceed tiniest target\n";
        return 6;
    }

    return 0;
}
