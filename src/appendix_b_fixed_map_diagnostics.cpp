#include "appendix_b_fixed_map_budget.hpp"
#include "appendix_b_multiplier_audit.hpp"
#include "appendix_b_parameter_scale.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char**argv){
    using namespace nsblowup;
    const std::string path=argc>1?argv[1]:"appendix_b_fixed_map.csv";
    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.04,161);
    const auto ls=appendix_b_lambda_scale_audit(p,0.05,4.1,0.05,20.0,0.05,0.04,161);
    const double h=std::exp(p.log_h);
    const double Lambda0=std::exp(ls.log_Lambda_for_U_tolerance);
    const double u0=std::exp(ls.max_log_abs_u0_Ymax);

    std::ofstream out(path);
    if(!out){std::cerr<<"failed to open output\n";return 1;}
    out<<std::setprecision(17);
    out<<"rho,lambda_factor,log_Lambda,chi_norm,zeta_norm,T_bound,inverse_one_plus_T_bound,"
          "phi_Wstar_Phi,phi_h_background_Phi,phi_B_etaAX_u_Phi,phi_B_detaAX_u_Phi,"
          "phi_h_u_Phi,phi_zeta_u_Phi,phi_Wstar_DXPhi,phi_B_etaAX_u_DXPhi,"
          "phi_B_detaAX_u_DXPhi,phi_Hstar_detaPhi,phi_u_detaPhi,"
          "u_linear,u_quadratic,u_Wstar_DXu,u_B_etaAX_u_DXu,u_B_detaAX_u_DXu,"
          "u_Hstar_detau,u_u_detau,u_pressure_p,u_pressure_peta,u_pressure_DXp,"
          "M_phi_map,M_u_map,contraction_bound,inverse_certified,contraction_certified\n";

    for(double rho:{1e-5,2e-5,5e-5,1e-4,2e-4}){
        const auto op=appendix_b_brho_operator_audit(18,6,rho);
        const auto mult=appendix_b_sample_multiplier_norms(h,0.05,sep.sigma_star,sep.eta0,rho,6,2401);
        for(double factor:{1.0,10.0,100.0,1000.0,10000.0}){
            const double Lambda=Lambda0*factor;
            const auto b=appendix_b_fixed_map_lipschitz_budget(
                op,mult.norms,h,Lambda,2.0,1.1*u0,1.0);
            out<<rho<<','<<factor<<','<<std::log(Lambda)<<','<<mult.norms.chi<<','<<mult.norms.zeta<<','
               <<b.T_bound<<','<<b.inverse_one_plus_T_bound<<','
               <<b.phi_Wstar_Phi<<','<<b.phi_h_background_Phi<<','
               <<b.phi_B_etaAX_u_Phi<<','<<b.phi_B_detaAX_u_Phi<<','
               <<b.phi_h_u_Phi<<','<<b.phi_zeta_u_Phi<<','
               <<b.phi_Wstar_DXPhi<<','<<b.phi_B_etaAX_u_DXPhi<<','
               <<b.phi_B_detaAX_u_DXPhi<<','<<b.phi_Hstar_detaPhi<<','<<b.phi_u_detaPhi<<','
               <<b.u_linear<<','<<b.u_quadratic<<','<<b.u_Wstar_DXu<<','
               <<b.u_B_etaAX_u_DXu<<','<<b.u_B_detaAX_u_DXu<<','
               <<b.u_Hstar_detau<<','<<b.u_u_detau<<','
               <<b.u_pressure_p<<','<<b.u_pressure_peta<<','<<b.u_pressure_DXp<<','
               <<b.M_phi_map<<','<<b.M_u_map<<','<<b.contraction_bound<<','
               <<(b.inverse_certified?1:0)<<','<<(b.contraction_certified?1:0)<<'\n';
        }
    }
    return 0;
}
