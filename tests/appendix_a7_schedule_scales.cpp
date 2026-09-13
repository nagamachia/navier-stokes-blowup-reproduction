#include "appendix_a7_schedule_scales.hpp"

#include <cmath>
#include <iostream>

int main() {
    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;

    const auto s = nsblowup::appendix_a7_schedule_scales(p, 0.4, 20.0, 0.05, 0.05);
    if (!(std::isfinite(s.log_X_patch) && std::isfinite(s.log_e_patch) &&
          std::isfinite(s.log_X_tail) && std::isfinite(s.log_e_tail))) {
        std::cerr << "non-finite A.7 schedule log scales\n";
        return 1;
    }
    if (!(s.log_Xpatch_over_Xtail < -1.0e6)) {
        std::cerr << "A.7 patch/tail radial separation is unexpectedly small\n";
        return 2;
    }
    if (!(s.log_epatch_over_etail > 5.0e5)) {
        std::cerr << "A.7 patch/tail amplitude separation is unexpectedly small\n";
        return 3;
    }
    if (!(s.log_target_scale_Cp < -7.0e2 && s.log_target_scale_I > 7.0e2)) {
        std::cerr << "A.7 target conversion should require log-space arithmetic\n";
        return 4;
    }

    return 0;
}
