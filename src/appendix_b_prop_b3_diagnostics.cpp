#include "appendix_b_analytic_multiplier.hpp"
#include "appendix_b_axis.hpp"
#include "appendix_b_axis_asymptotic.hpp"
#include "appendix_b_brho.hpp"
#include "appendix_b_parameter_scale.hpp"
#include "appendix_b_pressure_analytic.hpp"
#include "appendix_b_proposition_b2_certificate.hpp"
#include "appendix_b_proposition_b3_certificate.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char**argv){
    using namespace nsblowup;
    const std::string path=argc>1?argv[1]:"appendix_b_prop_b3_certificate.csv";
    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    constexpr double rho=1e-5;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.05,121);
    const auto scale=appendix_b_lambda_scale_audit(p,0.05,4.1,0.05,20.0,0.05,0.05,81);
    const auto center=appendix_b_Z_over_L_analytic_audit(p,0.05,rho,4*rho,8*rho,4.1,20.0,0.05,0.02);
    const auto op=appendix_b_brho_infinite_alpha_beta_operator_audit(rho);
    const double h=std::exp(p.log_h);
    const double Lambda0=std::exp(scale.log_Lambda_for_U_tolerance);
    const double u0=std::exp(center.log_u0_brho_Ymax);

    std::ofstream out(path);
    if(!out){ std::cerr<<"failed to open output\n"; return 1; }
    out<<std::setprecision(17);
    out<<"rho,lambda_factor,log_Lambda,b2_certificate,ball_radius,S0_Ymax,S0_exit,S1_exit,"
          "phi0_lower,phi_correction_sup,phi_lower,endpoint_phi_error,endpoint_DXphi_error,"
          "unperturbed_azimuthal_shear,perturbed_azimuthal_shear,endpoint_cex,"
          "log_C_chosen,log_C_required_axial,log_C_axial_margin,phi_positive,"
          "azimuthal_branch_certified,axial_branch_certified,B19_certified,proposition_b3_certificate\n";

    for(double factor=1.0;factor<=1e60;factor*=10.0){
        const double Lambda=Lambda0*factor;
        if(!std::isfinite(Lambda)) break;
        const auto am=appendix_b_analytic_multiplier_norms(h,0.05,sep.sigma_star,rho,0.45,Lambda);
        const auto b2=appendix_b_proposition_b2_certificate(op,am.norms,h,Lambda,u0,am.scalar_cauchy_factor,center.certified);
        AppendixBPropositionB3Certificate b3{};
        bool evaluated=false;
        if(b2.proposition_b2_certificate &&
           b2.ball_radius*appendix_b_brho_point_eval_S0(4.1)<0.265){
            const auto ep=appendix_b_endpoint_asymptotic_audit(p,sep,Lambda,0.05,241);
            b3=appendix_b_proposition_b3_certificate(b2,ep,am.log_C_for_g);
            evaluated=true;
        }
        out<<rho<<','<<factor<<','<<std::log(Lambda)<<','<<(b2.proposition_b2_certificate?1:0)<<','
           <<b2.ball_radius<<','<<appendix_b_brho_point_eval_S0(4.1)<<','
           <<appendix_b_brho_point_eval_S0(4.0)<<','<<appendix_b_brho_point_eval_S1(4.0)<<',';
        if(evaluated){
            out<<b3.phi0_global_lower<<','<<b3.phi_correction_sup<<','<<b3.phi_lower<<','
               <<b3.endpoint_phi_correction_sup<<','<<b3.endpoint_DXphi_correction_sup<<','
               <<b3.unperturbed_azimuthal_shear_lower<<','<<b3.perturbed_azimuthal_shear_lower<<','
               <<b3.endpoint_cex<<','<<b3.log_C_chosen<<','<<b3.log_C_required_axial<<','
               <<b3.log_C_axial_margin<<','<<(b3.phi_positive?1:0)<<','
               <<(b3.azimuthal_branch_certified?1:0)<<','<<(b3.axial_branch_certified?1:0)<<','
               <<(b3.endpoint_B19_certified?1:0)<<','<<(b3.proposition_b3_certificate?1:0);
        }else{
            out<<"nan,nan,nan,nan,nan,nan,nan,0,nan,nan,nan,0,0,0,0,0";
        }
        out<<'\n';
    }
    return 0;
}
