#pragma once

#include "appendix_a7_schedule_scales.hpp"

#include <algorithm>
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

namespace detail {

template<class F>
inline double positive_simpson_log_integral(F&& f, double y_max,
                                            std::size_t panels = 4096) {
    if (!(y_max > 0.5)) throw std::invalid_argument("A.7 log integral window is too short");
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    const double dy = y_max / static_cast<double>(panels);
    double sum = 0.0;
    for (std::size_t i = 0; i <= panels; ++i) {
        const double y = dy * static_cast<double>(i);
        const double w = (i == 0 || i == panels) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum += w * f(y);
    }
    const double value = sum * dy / 3.0;
    if (!(value > 0.0 && std::isfinite(value)))
        throw std::runtime_error("A.7 leading heat integral is not positive finite");
    return std::log(value);
}

} // namespace detail

// Leading large-X_K form obtained from
// H(2d/X)=1-2h(1+h)d/X+O(X^-2), d=1-eta^2.
// Keeping the result in signed-log form avoids underflow in the ordered regime.
inline AppendixA7LeadingHeatLogDiscrepancy
appendix_a7_leading_heat_log_discrepancy(
    const OuterProfileParameters& p, double eta,
    double Tf = 20.0, double co = 0.05,
    double y_max = 160.0, std::size_t panels = 4096) {
    const auto scales = appendix_a7_schedule_scales(p, eta, Tf, co);
    const double h = std::exp(p.log_h);
    const double d = 1.0 - eta * eta;
    if (!(h > 0.0 && d > 0.0))
        throw std::invalid_argument("leading A.7 heat log discrepancy requires h>0 and |eta|<1");

    const double log_a = std::log(2.0) + std::log(h) + std::log1p(h) + std::log(d);
    auto chi = [](double y) { return appendix_a_smooth_step((y - 0.2) / 0.3); };

    const double log_int_cp = detail::positive_simpson_log_integral(
        [&](double y) { return chi(y) * std::exp(-(2.0 + 2.0 * h) * y); },
        y_max, panels);
    const double log_int_s = detail::positive_simpson_log_integral(
        [&](double y) { return chi(y) * std::exp(-(1.0 + 2.0 * h) * y); },
        y_max, panels);
    const double log_int_i = detail::positive_simpson_log_integral(
        [&](double y) { return chi(y) * std::exp(-h * y); },
        y_max, panels);

    AppendixA7LeadingHeatLogDiscrepancy out;
    out.Cp = {-1, log_a + log_int_cp - scales.log_X_tail};
    out.S  = {+1, log_a + log_int_s  - scales.log_X_tail};
    out.I  = {-1, 0.5 * std::log(2.0) + log_a + log_int_i - scales.log_X_tail};

    // The compensation target has the opposite sign and is converted from
    // tail natural units to the second-patch normalized units.
    out.patch_target_Cp = {+1, out.Cp.log_abs + scales.log_target_scale_Cp};
    out.patch_target_S  = {-1, out.S.log_abs  + scales.log_target_scale_S};
    out.patch_target_I  = {+1, out.I.log_abs  + scales.log_target_scale_I};
    return out;
}

} // namespace nsblowup
