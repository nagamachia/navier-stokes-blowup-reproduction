#include "appendix_a3_global.hpp"

#include <cmath>
#include <iostream>

namespace {
bool is_finite_value(double x){return std::isfinite(x);}
}

int main(){
    using namespace nsblowup;

    OuterProfileParameters p;
    p.Md=2.0;
    p.log_Pstar=25.0;
    p.lambda=0.05;
    p.log_h=-40.0;
    p.validate();

    const auto a=appendix_a3_integrate_to_pulse(p,0.4,0.02);
    if (!(is_finite_value(a.m_at_pulse)&&is_finite_value(a.j_at_pulse)&&is_finite_value(a.s_at_pulse)&&
          is_finite_value(a.rI_at_pulse)&&is_finite_value(a.log_E_at_pulse))) {
        std::cerr<<"non-finite A.14 normalized state\n";
        return 1;
    }
    if (!(std::abs(a.m_at_pulse)<1e-3 && std::abs(a.j_at_pulse)<1e-3)) {
        std::cerr<<"pre-pulse M/J ratios did not decay: "<<a.m_at_pulse<<" "<<a.j_at_pulse<<"\n";
        return 1;
    }
    if (!(std::abs(a.pre_M_at_first_bump)<std::abs(a.m_at_pulse) &&
          std::abs(a.pre_J_at_first_bump)<std::abs(a.j_at_pulse))) {
        std::cerr<<"A.15 bump-center normalization did not further suppress discrepancies\n";
        return 1;
    }

    const double amp=appendix_a3_solve_amplitude(0.4);
    const auto closed=appendix_a3_close_MJ_from_schedule(p,0.4,amp,0.02);
    if (!(std::abs(closed.residual_M)<1e-12 && std::abs(closed.residual_J)<1e-12)) {
        std::cerr<<"schedule-driven M/J closure failed\n";
        return 1;
    }

    const auto b=appendix_a3_integrate_to_pulse(p,-0.4,0.02);
    if (!(std::abs(a.m_at_pulse+b.m_at_pulse)<1e-10 &&
          std::abs(a.j_at_pulse+b.j_at_pulse)<1e-10 &&
          std::abs(a.s_at_pulse-b.s_at_pulse)<1e-10 &&
          std::abs(a.rI_at_pulse-b.rI_at_pulse)<1e-12 &&
          std::abs(a.log_E_at_pulse-b.log_E_at_pulse)<1e-12)) {
        std::cerr<<"eta reflection symmetry failed\n";
        return 1;
    }

    OuterProfileParameters q=p;
    q.lambda=0.02;
    q.log_h=-40.0;
    q.validate();
    const auto c=appendix_a3_integrate_to_pulse(q,0.4,0.03);
    if (!(std::abs(c.m_at_pulse)<std::abs(a.m_at_pulse))) {
        std::cerr<<"lambda scaling of pre-pulse M ratio is inconsistent\n";
        return 1;
    }

    const auto low=appendix_a3_evaluate_global_closure(p,0.0,0.9,20.0,0.05,0.05);
    const auto high=appendix_a3_evaluate_global_closure(p,0.0,1.2,20.0,0.05,0.05);
    if (!(is_finite_value(low.s_infinity)&&is_finite_value(high.s_infinity)&&
          low.s_infinity<0.0&&high.s_infinity>0.0)) {
        std::cerr<<"global S(infinity) root is not bracketed: "<<low.s_infinity<<" "<<high.s_infinity<<"\n";
        return 1;
    }
    const auto global=appendix_a3_solve_global_amplitude(p,0.0,20.0,0.05,0.05,45);
    if (!(global.amplitude>0.9&&global.amplitude<1.2&&
          std::abs(global.s_infinity)<1e-8&&
          std::abs(global.pulse_residual_M)<1e-12&&
          std::abs(global.pulse_residual_J)<1e-12&&
          std::abs(global.rI_exterior_start-global.rI_target)<2e-5&&
          global.exterior_hold_length>0.0&&global.terminal_weight>0.0)) {
        std::cerr<<"global A.3 closure failed: amp="<<global.amplitude
                 <<" S="<<global.s_infinity
                 <<" rIerr="<<(global.rI_exterior_start-global.rI_target)<<"\n";
        return 1;
    }

    std::cout<<"Appendix A.2 -> A.3 global closure checks passed\n";
    return 0;
}
