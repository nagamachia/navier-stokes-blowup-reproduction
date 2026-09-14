#include "appendix_b_analytic_multiplier.hpp"
#include "appendix_b_axis.hpp"
#include "appendix_b_brho.hpp"
#include "appendix_b_parameter_scale.hpp"
#include "appendix_b_pressure_analytic.hpp"
#include "appendix_b_proposition_b2_certificate.hpp"

#include <cmath>
#include <iostream>

int main(){
    using namespace nsblowup;
    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    constexpr double rho=1e-5;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.05,121);
    const auto scale=appendix_b_lambda_scale_audit(p,0.05,4.1,0.05,20.0,0.05,0.05,81);
    const double h=std::exp(p.log_h);
    const double Lambda0=std::exp(scale.log_Lambda_for_U_tolerance);
    const auto op=appendix_b_brho_infinite_alpha_beta_operator_audit(rho);
    if(!op.infinite_eta_certified || !op.infinite_radial_certified) return 1;

    const auto center=appendix_b_Z_over_L_analytic_audit(
        p,0.05,rho,4.0*rho,8.0*rho,4.1,20.0,0.05,0.02);
    if(!center.certified || !center.pressure.certified) return 5;
    const double u0=std::exp(center.log_u0_brho_Ymax);
    if(!(u0>0.0) || !std::isfinite(u0)) return 6;

    const double sampled_u0=std::exp(scale.max_log_abs_u0_Ymax);
    if(!(u0>=sampled_u0)) return 7;

    AppendixBPropositionB2Certificate winner{};
    double winner_factor=0.0;
    for(double factor=1.0; factor<=1e60; factor*=10.0){
        const double Lambda=Lambda0*factor;
        if(!std::isfinite(Lambda)) break;
        const auto am=appendix_b_analytic_multiplier_norms(h,0.05,sep.sigma_star,rho,0.45,Lambda);
        const auto cert=appendix_b_proposition_b2_certificate(
            op,am.norms,h,Lambda,u0,am.scalar_cauchy_factor,center.certified);
        if(cert.proposition_b2_certificate){ winner=cert; winner_factor=factor; break; }
    }
    if(!(winner_factor>0.0)) return 2;
    if(!winner.infinite_alpha_beta || !winner.inverse_certified ||
       !winner.contraction_certified || !winner.invariant_ball_certified ||
       !winner.analytic_center_certified) return 3;
    if(!(winner.lipschitz.contraction_bound<1.0) ||
       !(winner.invariant_lhs<=winner.ball_radius) ||
       !(winner.invariant_margin>=0.0) ||
       !std::isfinite(winner.source.source_norm)) return 4;
    std::cout << "Appendix B.2 analytic-center certificate: factor=" << winner_factor
              << " log_u0=" << center.log_u0_brho_Ymax
              << " pi=" << center.pressure.pi0_brho_norm
              << " dpi=" << center.pressure.dpi0_brho_norm
              << " q=" << winner.lipschitz.contraction_bound
              << " source=" << winner.source.source_norm
              << " radius=" << winner.ball_radius
              << " lhs=" << winner.invariant_lhs
              << " margin=" << winner.invariant_margin << '\n';
    return 0;
}
