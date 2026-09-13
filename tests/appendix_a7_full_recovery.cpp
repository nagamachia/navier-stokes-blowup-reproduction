#include "appendix_a7_full_recovery.hpp"

#include <cmath>
#include <iostream>

int main() {
    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;

    for (double eta : {-0.8, -0.4, 0.0, 0.4, 0.8}) {
        const auto a = nsblowup::appendix_a7_full_recovery_audit(
            p, eta, 20.0, 0.05, 1024);

        if (!a.heat_asymptotic_certified) {
            std::cerr << "A.7 heat asymptotic bound not certified\n";
            return 1;
        }
        if (!(a.worst_heat_log_relative_bound < -1e6)) {
            std::cerr << "A.7 ordered heat error bound is unexpectedly weak\n";
            return 2;
        }
        if (!a.dominant_compensation_certified) {
            std::cerr << "A.7 dominant compensation hierarchy not certified\n";
            return 3;
        }
        if (!a.pressure_tracks_Cp) {
            std::cerr << "A.7 Cp/pressure coupling lost\n";
            return 4;
        }
        if (a.direct_double_exact_recovery_available) {
            std::cerr << "ordered A.7 regime must not be mislabeled as direct double exact recovery\n";
            return 5;
        }
        if (!(a.compensation.target_log_spread > 1e5)) {
            std::cerr << "A.7 ordered scale hierarchy unexpectedly collapsed\n";
            return 6;
        }
        if (!(a.leading.patch_target_Cp.sign == 1 &&
              a.leading.patch_target_S.sign == -1 &&
              a.leading.patch_target_I.sign == 1)) {
            std::cerr << "A.7 compensation signs changed\n";
            return 7;
        }
    }

    return 0;
}
