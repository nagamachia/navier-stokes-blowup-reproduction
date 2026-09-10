#pragma once

#include "appendix_a3_global.hpp"
#include "outer_profile_stages.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

// Appendix A.4, (A.21)--(A.23).
// We use y=log(X/X_R) and normalize pressure by P_*^2.  This avoids
// unnecessary huge magnitudes while preserving the exact scheduled datum.
// The A.11 angular bumps are intentionally omitted here: their construction
// preserves the total pressure increment, while the full A.11 interval is kept.
struct A4PressureDatum {
    double eta{};
    double normalized_pi0{};       // Pi0 / P_*^2
    double normalized_deta{};      // (partial_eta Pi0) / P_*^2
    double normalized_detaeta{};   // (partial_eta^2 Pi0) / P_*^2
    double inner_contribution{};   // -5/2 f(eta)^2, normalized by P_*^2
    double exterior_hold_length{};
    double terminal_Qp{};
};

namespace detail {

struct A4RadialState {
    double log_g{};   // log(E^2/P_*^2)
    double integral{}; // integral E^2/P_*^2 dy accumulated forward
};

inline double a4_exp_safe(double x) {
    if (x < std::log(std::numeric_limits<double>::min())) return 0.0;
    if (x > std::log(std::numeric_limits<double>::max()))
        throw std::overflow_error("Appendix A.4 normalized E^2 overflow");
    return std::exp(x);
}

inline A4RadialState a4_advance_constant_l(A4RadialState q,
                                            double l,
                                            double length) {
    if (length <= 0.0) return q;
    const double rate = 2.0 * l - 1.0; // d_y log(E^2)
    const double g0 = a4_exp_safe(q.log_g);
    if (g0 > 0.0) {
        if (std::abs(rate) < 1e-14)
            q.integral += g0 * length;
        else
            q.integral += g0 * std::expm1(rate * length) / rate;
    }
    q.log_g += rate * length;
    return q;
}

template<class LFunc>
A4RadialState a4_integrate_l(A4RadialState q,
                             double length,
                             LFunc&& l,
                             double max_step) {
    if (length <= 0.0) return q;
    const int steps = std::max(1, static_cast<int>(std::ceil(length / max_step)));
    const double ds = length / static_cast<double>(steps);
    auto rhs = [&](const A4RadialState& z, double y) {
        return A4RadialState{2.0 * l(y) - 1.0, a4_exp_safe(z.log_g)};
    };
    for (int n = 0; n < steps; ++n) {
        const double y = n * ds;
        const auto k1 = rhs(q, y);
        const A4RadialState q2{q.log_g + 0.5*ds*k1.log_g,
                               q.integral + 0.5*ds*k1.integral};
        const auto k2 = rhs(q2, y + 0.5*ds);
        const A4RadialState q3{q.log_g + 0.5*ds*k2.log_g,
                               q.integral + 0.5*ds*k2.integral};
        const auto k3 = rhs(q3, y + 0.5*ds);
        const A4RadialState q4{q.log_g + ds*k3.log_g,
                               q.integral + ds*k3.integral};
        const auto k4 = rhs(q4, y + ds);
        q.log_g += ds*(k1.log_g + 2*k2.log_g + 2*k3.log_g + k4.log_g)/6.0;
        q.integral += ds*(k1.integral + 2*k2.integral + 2*k3.integral + k4.integral)/6.0;
    }
    return q;
}

struct A4InterpolationIntegral {
    double value{};
    double deta{};
    double detaeta{};
    double log_g_end{};
};

inline A4InterpolationIntegral a4_integrate_a10(double log_g_start,
                                                 double eta,
                                                 double lambda,
                                                 double Tf,
                                                 int panels) {
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    const double dy = Tf / static_cast<double>(panels);
    const double J = std::log1p(eta*eta);
    const double Jeta = 2.0*eta/(1.0 + eta*eta);
    const double Jetaeta = 2.0*(1.0 - eta*eta)/std::pow(1.0 + eta*eta, 2);
    const double a = std::log(2.0) - J;

    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (int i = 0; i <= panels; ++i) {
        const double y = i * dy;
        const double theta = outer_theta_f(y, Tf);
        const double log_g = log_g_start + (-1.0 - 2.0*lambda)*y
                           + 2.0*(theta - 1.0)*a;
        const double g = a4_exp_safe(log_g);
        const double dlog = -2.0 * theta * Jeta;
        const double ratio2 = 4.0*theta*theta*Jeta*Jeta
                            - 2.0*theta*Jetaeta;
        const double w = (i == 0 || i == panels) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        s0 += w * g;
        s1 += w * dlog * g;
        s2 += w * ratio2 * g;
    }

    const double log_g_end = log_g_start + (-1.0 - 2.0*lambda)*Tf
                           - 2.0*a;
    return {s0*dy/3.0, s1*dy/3.0, s2*dy/3.0, log_g_end};
}

} // namespace detail

inline A4PressureDatum appendix_a4_pressure_datum(const OuterProfileParameters& p,
                                                   double eta,
                                                   double Tf = 20.0,
                                                   double co = 0.05,
                                                   double max_step = 0.02) {
    p.validate();
    if (!(eta >= -1.0 && eta <= 1.0))
        throw std::invalid_argument("eta must lie in [-1,1]");
    if (!(Tf > 0.0 && co > 0.0 && max_step > 0.0))
        throw std::invalid_argument("Tf, co and max_step must be positive");

    const double h = std::exp(p.log_h);
    if (!(h > 0.0))
        throw std::runtime_error("h underflowed in Appendix A.4 pressure evaluation");

    const double f = outer_shape_f(eta);
    const double Jeta = 2.0*eta/(1.0 + eta*eta);
    const double Jetaeta = 2.0*(1.0 - eta*eta)/std::pow(1.0 + eta*eta, 2);

    // y<=0: E/P_*=f exp(y/10), hence integral_{-inf}^0 E^2/P_*^2 dy=5 f^2.
    const double inner_integral = 5.0 * f * f;
    detail::A4RadialState q{2.0*std::log(f), inner_integral};

    // Section A.2 schedule up to the axial pulse.  E is unaffected by U.
    q = detail::a4_integrate_l(q, 1.0,
        [](double y) { return outer_l_first_transition(y); }, max_step);
    q = detail::a4_advance_constant_l(q, 0.0, p.Td());
    q = detail::a4_integrate_l(q, 1.0,
        [&](double y) { return outer_l_intermediate_entry(y, p.lambda); }, max_step);
    q = detail::a4_advance_constant_l(q, -p.lambda, p.Tw());
    q = detail::a4_advance_constant_l(q, -p.lambda, 13.0/p.lambda);

    // Before A.10 all eta dependence is exactly f(eta)^2.
    const double pre_a10_integral = q.integral;
    const double pre_deta = -2.0 * Jeta * pre_a10_integral;
    const double pre_detaeta = (4.0*Jeta*Jeta - 2.0*Jetaeta) * pre_a10_integral;

    // A.10 removes eta dependence continuously: E=c(y) f(eta)^theta(y).
    const int a10_panels = std::max(200, static_cast<int>(std::ceil(Tf/max_step/2.0))*2);
    const auto a10 = detail::a4_integrate_a10(q.log_g, eta, p.lambda, Tf, a10_panels);
    q.integral += a10.value;
    q.log_g = a10.log_g_end;

    // A.11 interval remains in the datum.  The two compact angular bumps are
    // omitted because the exact A.11 correction has zero total pressure increment.
    const double a11_length = 30.0 * std::log(1.0/p.lambda);
    q = detail::a4_advance_constant_l(q, -p.lambda, a11_length);

    // Exterior release and terminal collar, exactly as in the A.3 schedule.
    auto ldown = [&](double y) {
        return -p.lambda + (-1.0 + p.lambda) * appendix_a_smooth_step(y);
    };
    q = detail::a4_integrate_l(q, 1.0, ldown, max_step);

    double Q = (p.lambda - h)/(1.0 - p.lambda);
    Q = detail::integrate_Q(Q, 1.0, h, ldown, max_step);

    const double hold1 = 4.0 * std::log(1.0/h);
    q = detail::a4_advance_constant_l(q, -1.0, hold1);
    Q += (1.0 - h) * hold1;

    auto lup = [&](double y) {
        return -1.0 + (1.0 - h) * appendix_a_smooth_step(y);
    };
    q = detail::a4_integrate_l(q, 1.0, lup, max_step);
    Q = detail::integrate_Q(Q, 1.0, h, lup, max_step);

    const double rho_o = co * h;
    const double Qp = outer_terminal_Qp(h, rho_o, 600);
    if (!(Q > Qp && Qp > 0.0))
        throw std::runtime_error("Appendix A.4 exterior Q hold is not well ordered");
    const double hold_h = std::log(Q/Qp)/(1.0 - h);
    q = detail::a4_advance_constant_l(q, -h, hold_h);

    auto lterm = [&](double y) { return outer_terminal_l(y, h, rho_o); };
    q = detail::a4_integrate_l(q, 3.0, lterm, std::min(max_step, 0.01));

    // A.12 terminal power-law tail, y>=3.
    q.integral += detail::a4_exp_safe(q.log_g)/(1.0 + 2.0*h);

    A4PressureDatum out;
    out.eta = eta;
    out.normalized_pi0 = -0.5 * q.integral;
    out.normalized_deta = -0.5 * (pre_deta + a10.deta);
    out.normalized_detaeta = -0.5 * (pre_detaeta + a10.detaeta);
    out.inner_contribution = -2.5 * f * f;
    out.exterior_hold_length = hold_h;
    out.terminal_Qp = Qp;
    return out;
}

// (A.23) specialized to the exact inner branch y<=0.  It connects the datum
// at X=0 (y=-infinity) to the full pressure profile without introducing X_R:
// Pi(y)/P_*^2 = Pi0/P_*^2 + (5/2) f(eta)^2 exp(y/5).
inline double appendix_a4_inner_pressure_normalized(const A4PressureDatum& datum,
                                                     double y) {
    if (!(y <= 0.0))
        throw std::invalid_argument("inner pressure formula requires y<=0");
    const double f = outer_shape_f(datum.eta);
    return datum.normalized_pi0 + 2.5*f*f*std::exp(y/5.0);
}

inline double appendix_a4_scale_pressure(double normalized_value,
                                         double log_Pstar) {
    const double exponent = 2.0*log_Pstar;
    if (exponent > std::log(std::numeric_limits<double>::max()))
        throw std::overflow_error("P_*^2 pressure scale overflows double");
    return normalized_value * std::exp(exponent);
}

} // namespace nsblowup
