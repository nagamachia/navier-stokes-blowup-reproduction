#include "appendix_b_multiplier_audit.hpp"
#include "appendix_b_parameter_scale.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char**argv){
    using namespace nsblowup;
    const std::string path=argc>1?argv[1]:"appendix_b_brho.csv";
    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.04,161);
    const auto ls=appendix_b_lambda_scale_audit(p,0.05,4.1,0.05,20.0,0.05,0.04,161);
    const double h=std::exp(p.log_h), Lambda=std::exp(ls.log_Lambda_for_U_tolerance);
    const double u0=std::exp(ls.max_log_abs_u0_Ymax);

    std::ofstream out(path); if(!out){std::cerr<<"failed to open output\n";return 1;}
    out<<std::setprecision(17);
    out<<"rho,radial_order,eta_order,product_bound,DX_bound,I_bound,J1_bound,J2_bound,deta_bound,"
          "mixed_J1,mixed_J2,mixed_AX_J1,mixed_AX_J2,"
          "sigma_star,zeta_norm,invL_norm,Wstar_norm,Hstar_norm,"
          "M_r1_bracket_phi,M_r1_W_DXPhi,M_r1_H_detaPhi,"
          "M_r2_linear_u,M_r2_u2,M_r2_W_DXu,M_r2_H_deta_u,M_r2_u_deta_u,"
          "M_r2_pressure_p,M_r2_pressure_peta,M_r2_pressure_DXp,M_R1,M_R2,M,pressure_fraction,"
          "log_Lambda,beta_radial_proxy,log_contraction_proxy\n";
    // rho is chosen only after sigma_* in Appendix B.  The closest zeta_* poles
    // lie at a distance O(sigma_*/|H_*'|), about 5e-4 here, so sweep well below it.
    for(double rho:{1e-5,2e-5,5e-5,1e-4,2e-4,5e-4}){
        const auto op=appendix_b_brho_operator_audit(18,6,rho);
        const auto mult=appendix_b_sample_multiplier_norms(h,0.05,sep.sigma_star,sep.eta0,rho,6,2401);
        const auto b=appendix_b_remainder_lipschitz_budget(op,mult.norms,h,Lambda,2.0,1.1*u0,1.0);
        const double beta=std::max(op.J1/2.0,op.J2/2.0);
        const double logc=std::log(beta)+std::log(b.M)-std::log(Lambda);
        out<<rho<<','<<op.radial_order<<','<<op.eta_order<<','<<op.product<<','<<op.DX<<','<<op.I<<','<<op.J1<<','<<op.J2<<','<<op.deta<<','
           <<op.mixed_J1<<','<<op.mixed_J2<<','<<op.mixed_AX_J1<<','<<op.mixed_AX_J2<<','
           <<sep.sigma_star<<','<<mult.norms.zeta<<','<<mult.norms.invL<<','<<mult.norms.Wstar<<','<<mult.norms.Hstar<<','
           <<b.r1_bracket_phi<<','<<b.r1_W_DXPhi<<','<<b.r1_H_detaPhi<<','
           <<b.r2_linear_u<<','<<b.r2_u2<<','<<b.r2_W_DXu<<','<<b.r2_H_deta_u<<','<<b.r2_u_deta_u<<','
           <<b.r2_pressure_p<<','<<b.r2_pressure_peta<<','<<b.r2_pressure_DXp<<','<<b.M_R1<<','<<b.M_R2<<','<<b.M<<','<<b.pressure_fraction<<','
           <<ls.log_Lambda_for_U_tolerance<<','<<beta<<','<<logc<<'\n';
    }
    return 0;
}
