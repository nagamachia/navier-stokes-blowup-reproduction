#pragma once

#include "appendix_a3_global.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

// A.11 is a compact relative-E modification.  Its two coefficients are chosen
// from the exact local moment map; outside those supports the profile returns
// to the same l=-lambda power law.  Consequently it must update the cumulative
// moments, not be re-applied as a second global slope evolution.
inline double appendix_a3_a11_S_jump(double lambda,
                                     double c1,
                                     double c2,
                                     double weight_at_first_center,
                                     double span = 3.0,
                                     int panels = 2400) {
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    const double y1 = 30.0 * std::log(1.0 / lambda) - 3.0;
    const double y2 = y1 + 2.0;
    const double dy = span / static_cast<double>(panels);
    auto integrand = [&](double t) {
        const double y = y1 + t;
        const double factor = 1.0 +
            c1 * appendix_a3_unit_bump(y, y1) +
            c2 * appendix_a3_unit_bump(y, y2);
        // Difference from the unmodified E profile in dS/dy=-XE^2/2.
        const double base_weight = weight_at_first_center * std::exp(-2.0 * lambda * t);
        return -0.5 * base_weight * (factor * factor - 1.0);
    };
    double sum = integrand(0.0) + integrand(span);
    for (int i = 1; i < panels; ++i)
        sum += (i % 2 ? 4.0 : 2.0) * integrand(i * dy);
    return sum * dy / 3.0;
}

inline A3GlobalClosure appendix_a3_evaluate_global_closure_direct(
    const OuterProfileParameters& p,
    double eta,
    double amplitude,
    double Tf = 20.0,
    double co = 0.05,
    double max_step = 0.03) {
    p.validate();
    const double h = std::exp(p.log_h);
    if (!(h > 0.0))
        throw std::runtime_error("h underflowed; increase log_h for numerical evaluation");
    if (!(Tf > 0.0 && co > 0.0))
        throw std::invalid_argument("Tf and co must be positive");

    const auto pre = appendix_a3_integrate_to_pulse(p, eta, max_step);
    const auto mj = appendix_a3_close_pulse_MJ(
        p.lambda, pre.pre_M_at_first_bump, pre.pre_J_at_first_bump, amplitude);

    const double pulse_len = 13.0 / p.lambda;
    const double y1 = pulse_len - 3.0;
    const double y2 = pulse_len - 1.0;
    const int n = std::max(2, static_cast<int>(std::ceil(pulse_len / max_step / 2.0)) * 2);
    const double dy = pulse_len / static_cast<double>(n);
    auto pulse_integrand = [&](double y) {
        const double R = amplitude * appendix_a3_R0(p.lambda * y) +
                         mj.c1 * appendix_a3_unit_bump(y, y1) +
                         mj.c2 * appendix_a3_unit_bump(y, y2);
        return std::exp(-2.0 * p.lambda * y) * (R * R - 0.5);
    };
    double pulse_sum = pulse_integrand(0.0) + pulse_integrand(pulse_len);
    for (int i = 1; i < n; ++i)
        pulse_sum += (i % 2 ? 4.0 : 2.0) * pulse_integrand(i * dy);

    detail::SWeightState sw{
        pre.s_at_pulse + pulse_sum * dy / 3.0,
        std::exp(-26.0)};
    double rI = detail::advance_rI_constant_l(pre.rI_at_pulse, -p.lambda, pulse_len);

    // A.10: remove eta dependence while U remains zero.
    auto lint = [&](double y) {
        return detail::profile_interpolation_l(y, eta, p.lambda, Tf);
    };
    sw = detail::integrate_S_weight(sw, Tf, lint, max_step);
    A3NormalizedState rq{0, 0, 0, rI, 0};
    rq = detail::integrate_stage(rq, Tf, [](double) { return 0.0; }, lint, max_step);
    rI = rq.rI;

    // Constant A.10 profile up to the first A.11 correction support.
    const double uniform_len = 30.0 * std::log(1.0 / p.lambda);
    const double first_center = uniform_len - 3.0;
    rI = detail::advance_rI_constant_l(rI, -p.lambda, first_center);
    sw = detail::advance_S_constant_l(sw, -p.lambda, first_center);

    const double target_rI = 1.0 / (1.0 - p.lambda);
    // As noted in the proof, after conversion to the first-center
    // normalization the exact required I change is target_rI-rI.
    const auto ip = appendix_a3_close_angular_I_pressure(
        p.lambda, target_rI - rI, 0.0);

    // Apply A.11 as a compact cumulative-moment jump, then advance the
    // underlying l=-lambda profile across its three-unit support window.
    sw.sbar += appendix_a3_a11_S_jump(
        p.lambda, ip.c1, ip.c2, sw.weight, 3.0, 2400);
    sw = detail::advance_S_constant_l(sw, -p.lambda, 3.0);
    rI = target_rI;

    // Exterior release from -lambda to -1.
    auto ldown = [&](double y) {
        return -p.lambda + (-1.0 + p.lambda) * appendix_a_smooth_step(y);
    };
    sw = detail::integrate_S_weight(sw, 1.0, ldown, max_step);
    double Q = (p.lambda - h) / (1.0 - p.lambda);
    Q = detail::integrate_Q(Q, 1.0, h, ldown, max_step);

    // Exact l=-1 hold from A.17.
    const double hold1 = 4.0 * std::log(1.0 / h);
    sw = detail::advance_S_constant_l(sw, -1.0, hold1);
    Q += (1.0 - h) * hold1;

    auto lup = [&](double y) {
        return -1.0 + (1.0 - h) * appendix_a_smooth_step(y);
    };
    sw = detail::integrate_S_weight(sw, 1.0, lup, max_step);
    Q = detail::integrate_Q(Q, 1.0, h, lup, max_step);

    const double rho_o = co * h;
    const double Qp = outer_terminal_Qp(h, rho_o, 600);
    if (!(Q > Qp && Qp > 0.0))
        throw std::runtime_error("exterior Q hold is not well ordered");
    const double hold_h = std::log(Q / Qp) / (1.0 - h);
    sw = detail::advance_S_constant_l(sw, -h, hold_h);

    auto lterm = [&](double y) { return outer_terminal_l(y, h, rho_o); };
    sw = detail::integrate_S_weight(sw, 3.0, lterm, std::min(max_step, 0.01));

    // Exact infinite Epow tail.
    const double s_infinity = sw.sbar - sw.weight / (4.0 * h);

    A3GlobalClosure out;
    out.amplitude = amplitude;
    out.s_infinity = s_infinity;
    out.pulse_residual_M = mj.residual_M;
    out.pulse_residual_J = mj.residual_J;
    out.rI_exterior_start = rI;
    out.rI_target = target_rI;
    out.angular_c1 = ip.c1;
    out.angular_c2 = ip.c2;
    out.exterior_hold_length = hold_h;
    out.terminal_weight = sw.weight;
    return out;
}

inline A3GlobalClosure appendix_a3_solve_global_amplitude_direct(
    const OuterProfileParameters& p,
    double eta,
    double Tf = 20.0,
    double co = 0.05,
    double max_step = 0.03,
    int iterations = 48) {
    auto value = [&](double a) {
        return appendix_a3_evaluate_global_closure_direct(p, eta, a, Tf, co, max_step);
    };
    double lo = 0.9, hi = 1.2;
    auto vlo = value(lo), vhi = value(hi);
    if (!(vlo.s_infinity < 0.0 && vhi.s_infinity > 0.0))
        throw std::runtime_error("global A.3 S(infinity) root not bracketed on [0.9,1.2]");
    for (int i = 0; i < iterations; ++i) {
        const double mid = 0.5 * (lo + hi);
        const auto vm = value(mid);
        if (vm.s_infinity > 0.0) hi = mid;
        else lo = mid;
    }
    return value(0.5 * (lo + hi));
}

} // namespace nsblowup
