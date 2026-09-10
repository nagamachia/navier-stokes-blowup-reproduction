#include "appendix_a4_pressure.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace nsblowup;

namespace {

bool near(double a, double b, double atol, double rtol) {
    return std::abs(a-b) <= atol + rtol*std::max(std::abs(a),std::abs(b));
}

OuterProfileParameters theorem_scale_parameters() {
    OuterProfileParameters p;
    p.Md = 4.0;
    p.log_Pstar = 70.0;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.validate();
    return p;
}

} // namespace

int main() {
    const auto p = theorem_scale_parameters();
    constexpr double Tf = 20.0;
    constexpr double co = 0.05;
    constexpr double step = 0.03;

    // (A.21)/(A.22): Pi0 is even, contains the exact inner contribution
    // -(5/2)P_*^2 f^2, and is more negative after adding the exterior integral.
    for (double eta : {0.0, 0.25, 0.6, 0.9}) {
        const auto plus = appendix_a4_pressure_datum(p, eta, Tf, co, step);
        const auto minus = appendix_a4_pressure_datum(p, -eta, Tf, co, step);
        if (!std::isfinite(plus.normalized_pi0) ||
            !std::isfinite(plus.normalized_deta) ||
            !std::isfinite(plus.normalized_detaeta)) {
            std::cerr << "non-finite Appendix A.4 pressure datum at eta=" << eta << "\n";
            return 1;
        }
        if (!near(plus.normalized_pi0, minus.normalized_pi0, 2e-12, 2e-11)) {
            std::cerr << "Pi0 evenness failed at eta=" << eta << "\n";
            return 1;
        }
        if (!near(plus.normalized_deta, -minus.normalized_deta, 2e-12, 2e-10)) {
            std::cerr << "Pi0 derivative oddness failed at eta=" << eta << "\n";
            return 1;
        }
        if (!(plus.normalized_pi0 <= plus.inner_contribution + 1e-13)) {
            std::cerr << "A.22 inner pressure bound failed at eta=" << eta
                      << ": " << plus.normalized_pi0
                      << " > " << plus.inner_contribution << "\n";
            return 1;
        }
        if (eta > 0.0 && !(eta*plus.normalized_deta > 0.0)) {
            std::cerr << "A.22 pressure monotonicity failed at eta=" << eta << "\n";
            return 1;
        }
    }

    const auto axis = appendix_a4_pressure_datum(p, 0.0, Tf, co, step);
    if (std::abs(axis.normalized_deta) > 2e-12) {
        std::cerr << "Pi0_eta(0) must vanish by evenness\n";
        return 1;
    }

    // Differentiate the actual scheduled integral and compare with the analytic
    // eta derivatives used by Appendix B's Z_* datum.
    const double eta = 0.47;
    const double eps = 2e-5;
    const auto center = appendix_a4_pressure_datum(p, eta, Tf, co, step);
    const auto left = appendix_a4_pressure_datum(p, eta-eps, Tf, co, step);
    const auto right = appendix_a4_pressure_datum(p, eta+eps, Tf, co, step);
    const double fd1 = (right.normalized_pi0-left.normalized_pi0)/(2.0*eps);
    const double fd2 = (right.normalized_pi0-2.0*center.normalized_pi0+left.normalized_pi0)/(eps*eps);
    if (!near(fd1, center.normalized_deta, 2e-8, 2e-6)) {
        std::cerr << "Pi0_eta finite-difference mismatch: " << fd1
                  << " vs " << center.normalized_deta << "\n";
        return 1;
    }
    if (!near(fd2, center.normalized_detaeta, 2e-5, 2e-4)) {
        std::cerr << "Pi0_etaeta finite-difference mismatch: " << fd2
                  << " vs " << center.normalized_detaeta << "\n";
        return 1;
    }

    // (A.23), exact inner branch.  y->-infinity recovers Pi0 and pressure
    // increases toward the exterior because d_y Pi = E^2/2 > 0.
    const double p_far = appendix_a4_inner_pressure_normalized(center, -120.0);
    const double p_mid = appendix_a4_inner_pressure_normalized(center, -5.0);
    const double p_zero = appendix_a4_inner_pressure_normalized(center, 0.0);
    if (!near(p_far, center.normalized_pi0, 2e-10, 2e-10) ||
        !(center.normalized_pi0 < p_mid && p_mid < p_zero && p_zero < 0.0)) {
        std::cerr << "A.23 inner pressure profile regression failed\n";
        return 1;
    }

    // Scaling by P_*^2 is only a representation change.
    const double scaled = appendix_a4_scale_pressure(center.normalized_pi0, p.log_Pstar);
    if (!std::isfinite(scaled) || !(scaled < 0.0)) {
        std::cerr << "physical pressure scaling failed\n";
        return 1;
    }

    std::cout << "Appendix A.4 pressure datum regression passed\n"
              << "Pi0/P*^2(0)=" << axis.normalized_pi0 << "\n"
              << "Pi0/P*^2(" << eta << ")=" << center.normalized_pi0 << "\n"
              << "deta=" << center.normalized_deta
              << " detaeta=" << center.normalized_detaeta << "\n"
              << "terminal hold=" << center.exterior_hold_length << "\n";
    return 0;
}
