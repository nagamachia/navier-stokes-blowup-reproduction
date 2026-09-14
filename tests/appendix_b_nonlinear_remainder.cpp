#include "appendix_b_fixed_map_budget.hpp"
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
    for (const auto* q : {&out.B, &out.W, &out.Hc, &out.p, &out.p_eta, &out.R1, &out.R2}) {
        if (q->size() != n) return 1;
        for (double x : *q) if (!std::isfinite(x)) return 2;
    }
    if (std::abs(out.p[1] - in.g * in.g) > 1e-14) return 3;
    const double expected_peta1 = 2.0 * in.g * in.g_eta +
                                  2.0 * in.g * in.g * in.Phi_eta[0];
    if (std::abs(out.p_eta[1] - expected_peta1) > 1e-14) return 4;
    const auto ax = appendix_b_AX(in.u, n);
    if (std::abs(ax[1] - in.u[1] / 2.0) > 1e-15 ||
        std::abs(ax[2] - in.u[2] / 3.0) > 1e-15) return 5;
    const auto dx = appendix_b_DX(in.u, n);
    if (std::abs(dx[1] - in.u[1]) > 1e-15 ||
        std::abs(dx[2] - 2.0 * in.u[2]) > 1e-15) return 6;

    OuterProfileParameters p;
    p.lambda = 1e-5; p.log_h = -80.0; p.log_Pstar = 70.0;
    const auto scale = appendix_b_lambda_scale_audit(p, 0.05, 4.1, 0.05, 20.0, 0.05, 0.05, 81);
    if (!std::isfinite(scale.log_Lambda_for_U_tolerance) ||
        !(scale.log_Lambda_for_U_tolerance > 2.0 * p.log_Pstar)) return 7;

    constexpr double rho = 1e-5;
    const auto brho = appendix_b_brho_operator_audit(18, 6, rho);
    if (!(brho.product > 0.0 && brho.J1 > 0.0 && brho.J2 > 0.0 &&
          brho.mixed_AX_J1 > 0.0 && brho.mixed_AX_J2 > 0.0 &&
          brho.deta_J1 > 0.0 && brho.deta_J2 > 0.0)) return 8;
    const auto sep = appendix_b_separation_audit(p, 0.05, 20.0, 0.05, 0.05, 121);
    const auto mult = appendix_b_sample_multiplier_norms(
        std::exp(p.log_h), 0.05, sep.sigma_star, sep.eta0, rho, 6, 1201);
    if (!(mult.norms.invL > 0.0 && mult.norms.zeta > 0.0 && mult.norms.chi > 0.0) ||
        !std::isfinite(mult.norms.zeta) || !std::isfinite(mult.norms.chi)) return 9;

    const double Lambda = std::exp(scale.log_Lambda_for_U_tolerance);
    const double u0_bound = std::exp(scale.max_log_abs_u0_Ymax);
    const auto fixed = appendix_b_fixed_map_lipschitz_budget(
        brho, mult.norms, std::exp(p.log_h), Lambda, 2.0, 1.1*u0_bound, 1.0);
    if (!fixed.inverse_certified || !(fixed.inverse_one_plus_T_bound > 0.0) ||
        !(fixed.contraction_bound > 0.0) || !std::isfinite(fixed.contraction_bound)) return 10;
    if (!(fixed.inverse_finite_prefix_bound > 1.0 &&
          fixed.inverse_analytic_tail_bound >= 0.0 &&
          fixed.inverse_one_plus_T_bound >= fixed.inverse_finite_prefix_bound &&
          fixed.inverse_factorial_K > 0.0 &&
          fixed.inverse_envelope_ratio <= 1.0 + 1e-12)) return 11;
    if (!(fixed.phi_B_detaAX_u_DXPhi > 0.0 && fixed.u_B_detaAX_u_DXu > 0.0)) return 12;

    double closing_factor = 0.0;
    double closing_bound = 0.0;
    for (double factor = 1.0; factor <= 1e12; factor *= 10.0) {
        const auto trial = appendix_b_fixed_map_lipschitz_budget(
            brho, mult.norms, std::exp(p.log_h), factor*Lambda, 2.0, 1.1*u0_bound, 1.0);
        if (trial.contraction_certified) {
            closing_factor = factor;
            closing_bound = trial.contraction_bound;
            break;
        }
    }
    if (!(closing_factor > 0.0)) {
        std::cerr << "full fixed-map audit did not close through 1e12 Lambda factor\n";
        return 13;
    }

    std::cout << "Appendix B full fixed map: rho=" << rho
              << " crudeT=" << fixed.T_bound
              << " invFull=" << fixed.inverse_one_plus_T_bound
              << " invFinite=" << fixed.inverse_finite_prefix_bound
              << " invTail=" << fixed.inverse_analytic_tail_bound
              << " K=" << fixed.inverse_factorial_K
              << " envelopeRatio=" << fixed.inverse_envelope_ratio
              << " contraction=" << fixed.contraction_bound
              << " closingFactor=" << closing_factor
              << " closingBound=" << closing_bound << '\n';
    return 0;
}
