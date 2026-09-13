#include "appendix_b_nonlinear_remainder.hpp"
#include "appendix_b_parameter_scale.hpp"

#include <cmath>
#include <iostream>

int main() {
    using namespace nsblowup;
    constexpr std::size_t n = 8;

    AppendixBNonlinearRemainderInput in;
    in.eta = 0.2;
    in.h = 1e-6;
    in.Lambda = 1000.0;
    in.sigma_star = 0.01;
    in.g = 0.7;
    in.g_eta = -0.2;
    in.Phi = BPolynomial(n, 0.0);
    in.Phi_eta = BPolynomial(n, 0.0);
    in.u = BPolynomial(n, 0.0);
    in.u_eta = BPolynomial(n, 0.0);
    in.Phi[0] = 1.0;
    in.Phi[1] = -0.1;
    in.Phi_eta[0] = 0.03;
    in.u[1] = 0.4;
    in.u[2] = -0.05;
    in.u_eta[1] = -0.2;

    const auto out = appendix_b_nonlinear_remainder(in, n);
    for (const auto* p : {&out.B, &out.W, &out.Hc, &out.p, &out.p_eta, &out.R1, &out.R2}) {
        if (p->size() != n) return 1;
        for (double x : *p) if (!std::isfinite(x)) return 2;
    }

    if (std::abs(out.p[1] - in.g * in.g) > 1e-14) {
        std::cerr << "B.15 pressure primitive coefficient mismatch\n";
        return 3;
    }

    const double expected_peta1 = 2.0 * in.g * in.g_eta +
                                  2.0 * in.g * in.g * in.Phi_eta[0];
    if (std::abs(out.p_eta[1] - expected_peta1) > 1e-14) {
        std::cerr << "B.15 pressure eta derivative mismatch\n";
        return 4;
    }

    const auto ax = appendix_b_AX(in.u, n);
    if (std::abs(ax[1] - in.u[1] / 2.0) > 1e-15 ||
        std::abs(ax[2] - in.u[2] / 3.0) > 1e-15) return 5;

    const auto dx = appendix_b_DX(in.u, n);
    if (std::abs(dx[1] - in.u[1]) > 1e-15 ||
        std::abs(dx[2] - 2.0 * in.u[2]) > 1e-15) return 6;

    OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;
    const auto scale = appendix_b_lambda_scale_audit(p, 0.05, 4.1, 0.05, 20.0, 0.05, 0.05, 81);
    if (!std::isfinite(scale.log_Lambda_for_U_tolerance) ||
        !(scale.log_Lambda_for_U_tolerance > 2.0 * p.log_Pstar)) {
        std::cerr << "ordered Appendix B Lambda scale was not larger than P_*^2\n";
        return 7;
    }
    std::cout << "Appendix B ordered scale: logLambda(Ucorr<=0.05)="
              << scale.log_Lambda_for_U_tolerance
              << ", log(Lambda/P_*^2)=" << scale.log_Lambda_over_Pstar2
              << ", eta(max|Z*|)=" << scale.eta_at_max << '\n';

    return 0;
}
