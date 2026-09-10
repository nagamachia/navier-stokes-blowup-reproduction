#include "appendix_a3_global.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {
bool finite(double x){return std::isfinite(x);}
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
    if (!(finite(a.m_at_pulse)&&finite(a.j_at_pulse)&&finite(a.s_at_pulse)&&finite(a.log_E_at_pulse))) {
        std::cerr<<"non-finite A.14 normalized state\n";
        return 1;
    }
    // A.14 predicts strong decay of the normalized M and J discrepancies.
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

    // Reflection check: E is even in eta and U odd, hence M,J change sign
    // while S and log E are even.
    const auto b=appendix_a3_integrate_to_pulse(p,-0.4,0.02);
    if (!(std::abs(a.m_at_pulse+b.m_at_pulse)<1e-10 &&
          std::abs(a.j_at_pulse+b.j_at_pulse)<1e-10 &&
          std::abs(a.s_at_pulse-b.s_at_pulse)<1e-10 &&
          std::abs(a.log_E_at_pulse-b.log_E_at_pulse)<1e-12)) {
        std::cerr<<"eta reflection symmetry failed\n";
        return 1;
    }

    // Smaller lambda gives a longer A.9 hold and hence a smaller A.14 M ratio.
    OuterProfileParameters q=p;
    q.lambda=0.02;
    q.log_h=-40.0;
    q.validate();
    const auto c=appendix_a3_integrate_to_pulse(q,0.4,0.02);
    if (!(std::abs(c.m_at_pulse)<std::abs(a.m_at_pulse))) {
        std::cerr<<"lambda scaling of pre-pulse M ratio is inconsistent\n";
        return 1;
    }

    std::cout<<"Appendix A.2 -> A.3 global pre-pulse integration checks passed\n";
    return 0;
}
