#include "appendix_a3_global_closure.hpp"

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
    // Ordered asymptotic regime from (A.6): lambda is chosen only after the
    // fixed Md,Td,P* data.  The full A.19 regression therefore uses a genuinely
    // small lambda rather than the larger values appropriate to local tests.
    p.lambda=1.0e-5;
    p.log_h=-40.0;
    p.validate();

    const double step=0.05;
    const auto a=appendix_a3_integrate_to_pulse(p,0.4,step);
    if (!(is_finite_value(a.m_at_pulse)&&is_finite_value(a.j_at_pulse)&&is_finite_value(a.s_at_pulse)&&
          is_finite_value(a.rI_at_pulse)&&is_finite_value(a.log_E_at_pulse))) {
        std::cerr<<"non-finite A.14 normalized state\n";
        return 1;
    }
    if (!(std::abs(a.m_at_pulse)<1e-6 && std::abs(a.j_at_pulse)<1e-6)) {
        std::cerr<<"pre-pulse M/J ratios did not decay: "<<a.m_at_pulse<<" "<<a.j_at_pulse<<"\n";
        return 1;
    }
    if (!(std::abs(a.pre_M_at_first_bump)<=std::abs(a.m_at_pulse) &&
          std::abs(a.pre_J_at_first_bump)<=std::abs(a.j_at_pulse))) {
        std::cerr<<"A.15 bump-center normalization did not suppress discrepancies\n";
        return 1;
    }

    const double amp=appendix_a3_solve_amplitude(0.4);
    const auto closed=appendix_a3_close_MJ_from_schedule(p,0.4,amp,step);
    if (!(std::abs(closed.residual_M)<1e-12 && std::abs(closed.residual_J)<1e-12)) {
        std::cerr<<"schedule-driven M/J closure failed\n";
        return 1;
    }

    const auto b=appendix_a3_integrate_to_pulse(p,-0.4,step);
    if (!(std::abs(a.m_at_pulse+b.m_at_pulse)<1e-10 &&
          std::abs(a.j_at_pulse+b.j_at_pulse)<1e-10 &&
          std::abs(a.s_at_pulse-b.s_at_pulse)<1e-7 &&
          std::abs(a.rI_at_pulse-b.rI_at_pulse)<1e-12 &&
          std::abs(a.log_E_at_pulse-b.log_E_at_pulse)<1e-12)) {
        std::cerr<<"eta reflection symmetry failed\n";
        return 1;
    }

    const auto low=appendix_a3_evaluate_global_closure_direct(p,0.0,0.9,20.0,0.05,step);
    const auto high=appendix_a3_evaluate_global_closure_direct(p,0.0,1.2,20.0,0.05,step);
    if (!(is_finite_value(low.s_infinity)&&is_finite_value(high.s_infinity)&&
          low.s_infinity<0.0&&high.s_infinity>0.0)) {
        std::cerr<<"global S(infinity) root is not bracketed: "<<low.s_infinity<<" "<<high.s_infinity<<"\n";
        return 1;
    }

    // Compare the directly integrated endpoints against the principal A.19
    // expression after multiplication by lambda.  The remaining difference
    // should now be small on the scale lambda log(1/lambda).
    const double Kb=appendix_a3_Kb();
    const double c=(1.0-std::exp(-26.0))/4.0;
    const double err_low=p.lambda*low.s_infinity-(0.9*0.9*Kb-c);
    const double err_high=p.lambda*high.s_infinity-(1.2*1.2*Kb-c);
    if (!(std::abs(err_low)<0.02 && std::abs(err_high)<0.02)) {
        std::cerr<<"A.19 direct remainder too large: "<<err_low<<" "<<err_high<<"\n";
        return 1;
    }

    const auto global=appendix_a3_solve_global_amplitude_direct(p,0.0,20.0,0.05,step,34);
    if (!(global.amplitude>0.9&&global.amplitude<1.2&&
          std::abs(global.s_infinity)<1e-4&&
          std::abs(global.pulse_residual_M)<1e-12&&
          std::abs(global.pulse_residual_J)<1e-12&&
          std::abs(global.rI_exterior_start-global.rI_target)<1e-12&&
          global.exterior_hold_length>0.0&&global.terminal_weight>0.0)) {
        std::cerr<<"global A.3 closure failed: amp="<<global.amplitude
                 <<" S="<<global.s_infinity
                 <<" rIerr="<<(global.rI_exterior_start-global.rI_target)<<"\n";
        return 1;
    }

    std::cout<<"Appendix A.2 -> A.3 global closure checks passed\n";
    return 0;
}
