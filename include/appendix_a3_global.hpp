#pragma once

#include "appendix_a3_closure.hpp"
#include "outer_profile_stages.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct A3NormalizedState {
    // A.14-normalized cumulative quantities at the current radius.
    double m{};  // M/(X E)
    double j{};  // J/(X H E)
    double s{};  // S/(X E^2)
    double log_E{};
};

struct A3PrePulseData {
    double eta{};
    double lambda{};
    double m_at_pulse{};
    double j_at_pulse{};
    double s_at_pulse{};
    double pre_M_at_first_bump{};
    double pre_J_at_first_bump{};
    double log_E_at_pulse{};
};

namespace detail {

inline A3NormalizedState a3_rhs(const A3NormalizedState& q,
                                double U,
                                double l) {
    // d(log E)/dy=l-1/2.  Since H=sqrt(2X)E, direct differentiation gives
    // m' = U/E -(1/2+l)m,
    // j' = U/E -(1/2+2l)j,
    // s' = (U/E)^2 -1/2 -2ls.
    const double u_over_e = U == 0.0 ? 0.0 : U * std::exp(-q.log_E);
    return {
        u_over_e - (0.5 + l) * q.m,
        u_over_e - (0.5 + 2.0 * l) * q.j,
        u_over_e * u_over_e - 0.5 - 2.0 * l * q.s,
        l - 0.5};
}

inline A3NormalizedState add_scaled(const A3NormalizedState& a,
                                    const A3NormalizedState& b,
                                    double scale) {
    return {a.m + scale*b.m,
            a.j + scale*b.j,
            a.s + scale*b.s,
            a.log_E + scale*b.log_E};
}

template<class ULFunc, class LFunc>
A3NormalizedState integrate_stage(A3NormalizedState q,
                                  double length,
                                  ULFunc&& U,
                                  LFunc&& l,
                                  double max_step = 0.02) {
    if (length <= 0.0) return q;
    const int steps = std::max(1, static_cast<int>(std::ceil(length/max_step)));
    const double h = length/static_cast<double>(steps);
    for (int n=0; n<steps; ++n) {
        const double y = n*h;
        const double u1=U(y), l1=l(y);
        const auto k1=a3_rhs(q,u1,l1);

        const auto q2=add_scaled(q,k1,0.5*h);
        const double u2=U(y+0.5*h), l2=l(y+0.5*h);
        const auto k2=a3_rhs(q2,u2,l2);

        const auto q3=add_scaled(q,k2,0.5*h);
        const auto k3=a3_rhs(q3,u2,l2);

        const auto q4=add_scaled(q,k3,h);
        const double u4=U(y+h), l4=l(y+h);
        const auto k4=a3_rhs(q4,u4,l4);

        q.m += h*(k1.m+2*k2.m+2*k3.m+k4.m)/6.0;
        q.j += h*(k1.j+2*k2.j+2*k3.j+k4.j)/6.0;
        q.s += h*(k1.s+2*k2.s+2*k3.s+k4.s)/6.0;
        q.log_E += h*(k1.log_E+2*k2.log_E+2*k3.log_E+k4.log_E)/6.0;
    }
    return q;
}

} // namespace detail

// Integrate the complete A.2 schedule up to X_p, the start of the axial pulse.
// This evaluates the pre-pulse discrepancies in (A.14) directly from the
// staged profile rather than supplying them as manufactured inputs.
inline A3PrePulseData appendix_a3_integrate_to_pulse(const OuterProfileParameters& p,
                                                     double eta,
                                                     double max_step = 0.02) {
    p.validate();
    if (!(eta >= -1.0 && eta <= 1.0))
        throw std::invalid_argument("eta must lie in [-1,1]");

    const double f=outer_shape_f(eta);
    const double logE0=p.log_Pstar+std::log(f);
    const double invE=std::exp(-logE0);

    // Exact integrals of the temporary reference power law (A.7) on 0<X<=XR.
    A3NormalizedState q;
    q.m=4.0*eta*invE;
    q.j=2.5*eta*invE;
    q.s=16.0*eta*eta*invE*invE-5.0/12.0;
    q.log_E=logE0;

    // First one-unit transition after x=1.
    q=detail::integrate_stage(q,1.0,
        [=](double){return 4.0*eta;},
        [](double y){return outer_l_first_transition(y);},max_step);

    // Axial reduction stage.  Its local coordinate has length Td; k becomes
    // identically zero well before the endpoint as specified after (A.7).
    q=detail::integrate_stage(q,p.Td(),
        [&](double y){return outer_axial_U(y,eta,p.Md);},
        [](double){return 0.0;},max_step);

    // Entry into the intermediate power law and the Tw hold (A.9).
    q=detail::integrate_stage(q,1.0,
        [](double){return 0.0;},
        [&](double y){return outer_l_intermediate_entry(y,p.lambda);},max_step);
    q=detail::integrate_stage(q,p.Tw(),
        [](double){return 0.0;},
        [&](double){return -p.lambda;},max_step);

    const double y1=13.0/p.lambda-3.0;
    const double s1=0.5-p.lambda;
    const double s2=0.5-2.0*p.lambda;

    A3PrePulseData out;
    out.eta=eta;
    out.lambda=p.lambda;
    out.m_at_pulse=q.m;
    out.j_at_pulse=q.j;
    out.s_at_pulse=q.s;
    // Convert the A.14 normalization at Xp to the first end-bump
    // normalizations used in (A.15).
    out.pre_M_at_first_bump=q.m*std::exp(-s1*y1);
    out.pre_J_at_first_bump=q.j*std::exp(-s2*y1);
    out.log_E_at_pulse=q.log_E;
    return out;
}

inline PulseMomentClosure appendix_a3_close_MJ_from_schedule(
    const OuterProfileParameters& p,
    double eta,
    double amplitude,
    double max_step = 0.02) {
    const auto pre=appendix_a3_integrate_to_pulse(p,eta,max_step);
    return appendix_a3_close_pulse_MJ(p.lambda,
                                      pre.pre_M_at_first_bump,
                                      pre.pre_J_at_first_bump,
                                      amplitude);
}

} // namespace nsblowup
