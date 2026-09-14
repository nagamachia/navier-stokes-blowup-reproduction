#include "appendix_b_analytic_multiplier.hpp"
#include "appendix_b_axis.hpp"
#include "appendix_b_axis_asymptotic.hpp"
#include "appendix_b_brho.hpp"
#include "appendix_b_parameter_scale.hpp"
#include "appendix_b_pressure_analytic.hpp"
#include "appendix_b_proposition_b2_certificate.hpp"
#include "appendix_b_proposition_b3_certificate.hpp"

#include <cmath>
#include <iostream>

int main(){
    using namespace nsblowup;
    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    constexpr double rho=1e-5;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.05,121);
    const auto scale=appendix_b_lambda_scale_audit(p,0.05,4.1,0.05,20.0,0.05,0.05,81);
    const auto center=appendix_b_Z_over_L_analytic_audit(p,0.05,rho,4*rho,8*rho,4.1,20.0,0.05,0.02);
    const auto op=appendix_b_brho_infinite_alpha_beta_operator_audit(rho);
    const double h=std::exp(p.log_h);
    const double Lambda0=std::exp(scale.log_Lambda_for_U_tolerance);
    const double u0=std::exp(center.log_u0_brho_Ymax);

    AppendixBPropositionB3Certificate winner{};
    double winner_factor=0.0;
    for(double factor=1.0;factor<=1e60;factor*=10.0){
        const double Lambda=Lambda0*factor;
        if(!std::isfinite(Lambda)) break;
        const auto am=appendix_b_analytic_multiplier_norms(h,0.05,sep.sigma_star,rho,0.45,Lambda);
        const auto b2=appendix_b_proposition_b2_certificate(op,am.norms,h,Lambda,u0,am.scalar_cauchy_factor,center.certified);
        if(!b2.proposition_b2_certificate) continue;
        // Positivity must first become small enough to inherit from f0>0.265.
        if(b2.ball_radius*appendix_b_brho_point_eval_S0(4.1)>=0.265) continue;
        const auto ep=appendix_b_endpoint_asymptotic_audit(p,sep,Lambda,0.05,241);
        const auto b3=appendix_b_proposition_b3_certificate(b2,ep,am.log_C_for_g);
        if(b3.proposition_b3_certificate){ winner=b3; winner_factor=factor; break; }
    }
    if(!(winner_factor>0.0)) return 1;
    if(!winner.inherited_from_b2 || !winner.phi_positive ||
       !winner.azimuthal_branch_certified || !winner.axial_branch_certified ||
       !winner.endpoint_B19_certified || !winner.proposition_b3_certificate) return 2;
    if(!(winner.phi_lower>0.0) || !(winner.perturbed_azimuthal_shear_lower>2.0) ||
       !(winner.endpoint_cex>0.0) || !(winner.log_C_axial_margin>0.0)) return 3;
    std::cout << "Appendix B.3 certificate: factor=" << winner_factor
              << " phi_lower=" << winner.phi_lower
              << " exit=" << winner.perturbed_azimuthal_shear_lower
              << " cex=" << winner.endpoint_cex
              << " logC_margin=" << winner.log_C_axial_margin << '\n';
    return 0;
}
