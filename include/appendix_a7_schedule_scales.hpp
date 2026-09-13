#pragma once

#include "appendix_a3_global.hpp"

#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct AppendixA7ScheduleScales {
    double log_X_patch{};
    double log_e_patch{};
    double log_X_tail{};
    double log_e_tail{};
    double log_Xpatch_over_Xtail{};
    double log_epatch_over_etail{};
    double log_target_scale_Cp{};
    double log_target_scale_S{};
    double log_target_scale_I{};
};

namespace detail {

template<class LFunc>
inline double integrate_log_e(double log_e, double length, LFunc&& l,
                              double max_step = 0.02) {
    if (length <= 0.0) return log_e;
    const int steps = std::max(1, static_cast<int>(std::ceil(length / max_step)));
    const double h = length / static_cast<double>(steps);
    for (int n = 0; n < steps; ++n) {
        const double y0 = n * h;
        const double y1 = y0 + 0.5 * h;
        const double y2 = y0 + h;
        log_e += h * ((l(y0) - 0.5) + 4.0 * (l(y1) - 0.5) +
                      (l(y2) - 0.5)) / 6.0;
    }
    return log_e;
}

} // namespace detail

// Expose the relative natural scales needed by Proposition A.7 without ever
// exponentiating the exponentially separated radial locations.  X is tracked
// only through log X.  On the heat patch the normalized map uses
// E=e_* f(eta) x^{-1/2-lambda}, while at the terminal tail E=e_K x^{-1/2-h}.
inline AppendixA7ScheduleScales appendix_a7_schedule_scales(
    const OuterProfileParameters& p, double eta,
    double Tf = 20.0, double co = 0.05, double max_step = 0.03) {
    p.validate();
    if (std::abs(eta) > 1.0) throw std::invalid_argument("eta must lie in [-1,1]");
    if (!(Tf > 0.0 && co > 0.0 && max_step > 0.0))
        throw std::invalid_argument("A.7 schedule scale parameters must be positive");

    const double h = std::exp(p.log_h);
    if (!(h > 0.0)) throw std::runtime_error("h underflowed in A.7 schedule scale audit");

    double logX = 0.0;
    double logE = p.log_Pstar + std::log(outer_shape_f(eta));

    auto advance_variable = [&](double length, auto&& l) {
        logE = detail::integrate_log_e(logE, length, l, max_step);
        logX += length;
    };
    auto advance_constant = [&](double length, double l) {
        logE += (l - 0.5) * length;
        logX += length;
    };

    advance_variable(1.0, [](double y) { return outer_l_first_transition(y); });
    advance_constant(p.Td(), 0.0);
    advance_variable(1.0, [&](double y) { return outer_l_intermediate_entry(y, p.lambda); });

    const auto patches = reserved_patches_a9(p.lambda);
    const double hold_start_X = logX;
    const double hold_start_E = logE;
    const double patch_offset = patches.heat_a;
    const double log_X_patch = hold_start_X + patch_offset;
    const double log_E_patch = hold_start_E + (-0.5 - p.lambda) * patch_offset;
    const double log_e_patch = log_E_patch - std::log(outer_shape_f(eta));

    advance_constant(p.Tw(), -p.lambda);

    const double pulse_len = 13.0 / p.lambda;
    advance_constant(pulse_len, -p.lambda);

    auto lint = [&](double y) {
        return detail::profile_interpolation_l(y, eta, p.lambda, Tf);
    };
    advance_variable(Tf, lint);

    const double uniform_len = 30.0 * std::log(1.0 / p.lambda);
    advance_constant(uniform_len, -p.lambda);

    auto ldown = [&](double y) {
        return -p.lambda + (-1.0 + p.lambda) * appendix_a_smooth_step(y);
    };
    advance_variable(1.0, ldown);

    const double hold1 = 4.0 * std::log(1.0 / h);
    advance_constant(hold1, -1.0);

    auto lup = [&](double y) {
        return -1.0 + (1.0 - h) * appendix_a_smooth_step(y);
    };
    advance_variable(1.0, lup);

    double Q = (p.lambda - h) / (1.0 - p.lambda);
    Q = detail::integrate_Q(Q, 1.0, h, ldown, max_step);
    Q += (1.0 - h) * hold1;
    Q = detail::integrate_Q(Q, 1.0, h, lup, max_step);
    const double rho_o = co * h;
    const double Qp = outer_terminal_Qp(h, rho_o, 600);
    if (!(Q > Qp && Qp > 0.0))
        throw std::runtime_error("exterior Q hold is not well ordered in A.7 scale audit");
    const double hold_h = std::log(Q / Qp) / (1.0 - h);
    advance_constant(hold_h, -h);

    auto lterm = [&](double y) { return outer_terminal_l(y, h, rho_o); };
    advance_variable(3.0, lterm);

    AppendixA7ScheduleScales out;
    out.log_X_patch = log_X_patch;
    out.log_e_patch = log_e_patch;
    out.log_X_tail = logX;
    out.log_e_tail = logE;
    out.log_Xpatch_over_Xtail = log_X_patch - logX;
    out.log_epatch_over_etail = log_e_patch - logE;

    // Multipliers taking a heat discrepancy expressed in tail natural units
    // to the target seen by the normalized second-patch solver.
    out.log_target_scale_Cp = -2.0 * out.log_epatch_over_etail;
    out.log_target_scale_S = -out.log_Xpatch_over_Xtail -
                             2.0 * out.log_epatch_over_etail;
    out.log_target_scale_I = -1.5 * out.log_Xpatch_over_Xtail -
                             out.log_epatch_over_etail;
    return out;
}

} // namespace nsblowup
