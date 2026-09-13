#pragma once

#include "appendix_b_axis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBLambdaScaleAudit {
    double max_log_abs_Zstar{};
    double eta_at_max{};
    double max_log_abs_u0_Ymax{};
    double log_Lambda_for_U_tolerance{};
    double log_Lambda_over_Pstar2{};
    double U_tolerance{};
    double Ymax{};
};

// B.13 has U=U*+Lambda^{-1}u with u0=-Y Z*/(2L).
// This computes the explicit leading-order Lambda scale needed to make
// max_{eta,0<=Y<=Ymax}|u0|/Lambda <= U_tolerance.
// It is a transparent scale requirement, not the manuscript's full B_rho
// contraction constant (whose generic constants are not numerically specified).
inline AppendixBLambdaScaleAudit appendix_b_lambda_scale_audit(
    const OuterProfileParameters& p,
    double U_tolerance = 0.05,
    double Ymax = 4.1,
    double j0 = 0.05,
    double Tf = 20.0,
    double co = 0.05,
    double max_step = 0.04,
    int eta_samples = 161) {
    p.validate();
    if (!(U_tolerance > 0.0) || !(Ymax > 0.0) || eta_samples < 3)
        throw std::invalid_argument("invalid Appendix B Lambda-scale audit parameters");

    const double h = std::exp(p.log_h);
    double max_log_Z = -std::numeric_limits<double>::infinity();
    double eta_at_max = 0.0;
    double max_log_u0 = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < eta_samples; ++i) {
        const double eta = -1.0 + 2.0 * static_cast<double>(i) /
                                      static_cast<double>(eta_samples - 1);
        const auto a = appendix_b_axis_data(p, eta, j0, Tf, co, max_step);
        const double az = std::abs(a.normalized_Zstar);
        if (!(az > 0.0) || !std::isfinite(az)) continue;
        const double log_Z = 2.0 * p.log_Pstar + std::log(az);
        const double L = appendix_b_L(eta, h);
        const double log_u0 = std::log(Ymax) + log_Z - std::log(2.0 * L);
        if (log_Z > max_log_Z) {
            max_log_Z = log_Z;
            eta_at_max = eta;
        }
        max_log_u0 = std::max(max_log_u0, log_u0);
    }

    if (!std::isfinite(max_log_u0))
        throw std::runtime_error("failed Appendix B Lambda-scale sweep");

    AppendixBLambdaScaleAudit out;
    out.max_log_abs_Zstar = max_log_Z;
    out.eta_at_max = eta_at_max;
    out.max_log_abs_u0_Ymax = max_log_u0;
    out.log_Lambda_for_U_tolerance = max_log_u0 - std::log(U_tolerance);
    out.log_Lambda_over_Pstar2 = out.log_Lambda_for_U_tolerance - 2.0 * p.log_Pstar;
    out.U_tolerance = U_tolerance;
    out.Ymax = Ymax;
    return out;
}

} // namespace nsblowup
