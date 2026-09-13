#include "appendix_b_multiplier_audit.hpp"
#include "appendix_b_nonlinear_remainder.hpp"
#include "appendix_b_parameter_scale.hpp"

#include <algorithm>
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

    // Finite (alpha,beta) audit of the manuscript B-rho weight (B.4).
    const auto brho = appendix_b_brho_operator_audit(18, 6, 0.01);
    if (!(brho.product > 0.0 && brho.J1 > 0.0 && brho.J2 > 0.0 && brho.deta > 0.0))
        return 8;
    const auto sep = appendix_b_separation_audit(p, 0.05, 20.0, 0.05, 0.05, 121);
    const auto mult = appendix_b_sample_multiplier_norms(
        std::exp(p.log_h), 0.05, sep.sigma_star, sep.eta0, 0.01, 6, 1201);
    if (!(mult.norms.invL > 0.0 && mult.norms.zeta > 0.0) ||
        !std::isfinite(mult.norms.zeta)) return 9;

    const double Lambda = std::exp(scale.log_Lambda_for_U_tolerance);
    const double u0_bound = std::exp(scale.max_log_abs_u0_Ymax);
    const auto budget = appendix_b_remainder_lipschitz_budget(
        brho, mult.norms, std::exp(p.log_h), Lambda, 2.0, 1.1*u0_bound, 1.0);
    if (!(budget.M > 0.0 && std::isfinite(budget.M))) return 10;
    const double beta_radial_proxy = std::max(brho.J1/2.0, brho.J2/2.0);
    const double log_contraction_proxy = std::log(beta_radial_proxy) +
        std::log(budget.M) - scale.log_Lambda_for_U_tolerance;
    if (!(log_contraction_proxy < 0.0)) {
        std::cerr << "finite B-rho remainder budget did not close\n";
        return 11;
    }

    std::cout << "Appendix B ordered scale: logLambda(Ucorr<=0.05)="
              << scale.log_Lambda_for_U_tolerance
              << ", log(Lambda/P_*^2)=" << scale.log_Lambda_over_Pstar2
              << ", eta(max|Z*|)=" << scale.eta_at_max
              << ", B-rho M=" << budget.M
              << ", log contraction proxy=" << log_contraction_proxy << '\n';

    return 0;
}
