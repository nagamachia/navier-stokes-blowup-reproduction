#pragma once

#include "appendix_a7_schedule_scales.hpp"

#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct AppendixA7SignedLogMoment {
    int sign{};
    double log_abs{};
};

struct AppendixA7LeadingHeatLogDiscrepancy {
    AppendixA7SignedLogMoment Cp{};
    AppendixA7SignedLogMoment S{};
    AppendixA7SignedLogMoment I{};
    AppendixA7SignedLogMoment patch_target_Cp{};
    AppendixA7SignedLogMoment patch_target_S{};
    AppendixA7SignedLogMoment patch_target_I{};
};

struct AppendixA7LeadingHeatErrorAudit {
    double log_relative_bound_CpS{};
    double log_relative_bound_I{};
};

namespace detail {

inline double cutoff_exp_integral_log(double alpha, std::size_t panels = 1024) {
    if (!(alpha > 0.0)) throw std::invalid_argument("A.7 decay rate must be positive");
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    constexpr double b = 0.5;
    const double dy = b / static_cast<double>(panels);
    double sum = 0.0;
    for (std::size_t i = 0; i <= panels; ++i) {
        const double y = dy * static_cast<double>(i);
        const double chi = appendix_a_smooth_step((y - 0.2) / 0.3);
        const double v = chi * std::exp(-alpha * y);
        const double w = (i == 0 || i == panels) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum += w * v;
    }
    const double transition = sum * dy / 3.0;
    const double tail = std::exp(-alpha * b) / alpha;
    return std::log(transition + tail);
}

} // namespace detail

// Leading large-X_K form obtained from
// H(2d/X)=1-2h(1+h)d/X+O(X^-2), d=1-eta^2.
// Keeping the result in signed-log form avoids underflow in the ordered regime.
inline AppendixA7LeadingHeatLogDiscrepancy
appendix_a7_leading_heat_log_discrepancy(
    const OuterProfileParameters& p, double eta,
    double Tf = 20.0, double co = 0.05,
    std::size_t panels = 1024) {
    const auto scales = appendix_a7_schedule_scales(p, eta, Tf, co);
    const double h = std::exp(p.log_h);
    const double d = 1.0 - eta * eta;
    if (!(h > 0.0 && d > 0.0))
        throw std::invalid_argument("leading A.7 heat log discrepancy requires h>0 and |eta|<1");

    const double log_a = std::log(2.0) + std::log(h) + std::log1p(h) + std::log(d);
    const double log_int_cp = detail::cutoff_exp_integral_log(2.0 + 2.0 * h, panels);
    const double log_int_s = detail::cutoff_exp_integral_log(1.0 + 2.0 * h, panels);
    const double log_int_i = detail::cutoff_exp_integral_log(h, panels);

    AppendixA7LeadingHeatLogDiscrepancy out;
    out.Cp = {-1, log_a + log_int_cp - scales.log_X_tail};
    out.S  = {+1, log_a + log_int_s  - scales.log_X_tail};
    out.I  = {-1, 0.5 * std::log(2.0) + log_a + log_int_i - scales.log_X_tail};

    out.patch_target_Cp = {+1, out.Cp.log_abs + scales.log_target_scale_Cp};
    out.patch_target_S  = {-1, out.S.log_abs  + scales.log_target_scale_S};
    out.patch_target_I  = {+1, out.I.log_abs  + scales.log_target_scale_I};
    return out;
}

// Certified leading-expansion audit.  From Taylor's theorem for the positive
// heat-factor integral,
//   0 <= H(Z)-[1-h(1+h)Z]
//      <= 1/2 h(1+h)^2(2+h) Z^2.
// Thus the relative error in the linear E discrepancy is bounded by
// (1+h)(2+h)d/X.  For Cp and S, the extra (delta E)^2 term adds at most
// h(1+h)d/X, giving the convenient uniform bound 2(1+h)^2 d/X.
// Since X >= X_tail throughout the replacement, evaluating at X_tail bounds
// the complete cutoff-weighted moment integrals as well.
inline AppendixA7LeadingHeatErrorAudit appendix_a7_leading_heat_error_audit(
    const OuterProfileParameters& p, double eta,
    double Tf = 20.0, double co = 0.05) {
    const auto scales = appendix_a7_schedule_scales(p, eta, Tf, co);
    const double h = std::exp(p.log_h);
    const double d = 1.0 - eta * eta;
    if (!(h > 0.0 && d > 0.0))
        throw std::invalid_argument("A.7 heat error audit requires h>0 and |eta|<1");

    AppendixA7LeadingHeatErrorAudit out;
    out.log_relative_bound_CpS = std::log(2.0) + 2.0 * std::log1p(h) +
                                 std::log(d) - scales.log_X_tail;
    out.log_relative_bound_I = std::log1p(h) + std::log(2.0 + h) +
                               std::log(d) - scales.log_X_tail;
    return out;
}

} // namespace nsblowup
