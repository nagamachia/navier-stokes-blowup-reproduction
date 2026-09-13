#include "appendix_b_multiplier_audit.hpp"
#include "appendix_b_parameter_scale.hpp"

#include <cmath>
#include <iostream>

int main(){
    using namespace nsblowup;
    const auto op=appendix_b_brho_operator_audit(18,6,0.01);
    if(!(op.product>0.0 && std::isfinite(op.product))) return 1;
    if(!(op.J1>0.0 && op.J2>0.0 && op.I>0.0)) return 2;

    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.04,161);
    const double h=std::exp(p.log_h);
    const auto mult=appendix_b_sample_multiplier_norms(h,0.05,sep.sigma_star,sep.eta0,0.01,6,1601);
    if(!(mult.norms.invL>0.0 && mult.norms.zeta>0.0) || !std::isfinite(mult.norms.zeta)) return 3;

    const auto ls=appendix_b_lambda_scale_audit(p,0.05,4.1,0.05,20.0,0.05,0.04,161);
    const double Lambda=std::exp(ls.log_Lambda_for_U_tolerance);
    const double u0=std::exp(ls.max_log_abs_u0_Ymax);
    const auto B=appendix_b_remainder_lipschitz_budget(op,mult.norms,h,Lambda,2.0,1.1*u0,1.0);
    if(!(B.M>0.0 && std::isfinite(B.M))) return 4;
    if(!(B.M_R1>0.0 && B.M_R2>0.0)) return 5;
    if(!(B.pressure_fraction>=0.0 && B.pressure_fraction<=1.0)) return 6;

    // The ordered Lambda scale is exponentially larger than every finite
    // sampled finite-truncation remainder budget in this regression.
    const double beta=std::max(op.J1/2.0,op.J2/2.0);
    const double contraction=beta*B.M/Lambda;
    if(!(contraction<1.0)) {
        std::cerr << "finite B-rho contraction audit did not close: " << contraction << '\n';
        return 7;
    }
    return 0;
}
