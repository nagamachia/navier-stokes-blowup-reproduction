#include "appendix_b_analytic_multiplier.hpp"
#include "appendix_b_axis.hpp"
#include "appendix_b_brho.hpp"
#include "appendix_b_parameter_scale.hpp"
#include "appendix_b_proposition_b2_certificate.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char**argv){
    using namespace nsblowup;
    const std::string path=argc>1?argv[1]:"appendix_b_prop_b2_certificate.csv";
    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    constexpr double rho=1e-5;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.05,121);
    const auto scale=appendix_b_lambda_scale_audit(p,0.05,4.1,0.05,20.0,0.05,0.05,81);
    const double h=std::exp(p.log_h);
    const double Lambda0=std::exp(scale.log_Lambda_for_U_tolerance);
    const double u0=std::exp(scale.max_log_abs_u0_Ymax);
    const auto op=appendix_b_brho_infinite_alpha_beta_operator_audit(rho);

    std::ofstream out(path);
    if(!out){ std::cerr<<"failed to open output\n"; return 1; }
    out<<std::setprecision(17);
    out<<"rho,lambda_factor,log_Lambda,phi0_bound,u0_bound,ball_radius,phi_domain,u_domain,"
          "source_phi,source_u,source_norm,lipschitz_q,invariant_lhs,invariant_margin,"
          "infinite_alpha_beta,inverse_certified,contraction_certified,invariant_ball_certified,proposition_b2_certificate\n";
    for(double factor=1.0;factor<=1e24;factor*=10.0){
        const double Lambda=Lambda0*factor;
        if(!std::isfinite(Lambda)) break;
        const auto am=appendix_b_analytic_multiplier_norms(h,0.05,sep.sigma_star,rho,0.45,Lambda);
        const auto c=appendix_b_proposition_b2_certificate(op,am.norms,h,Lambda,u0,am.scalar_cauchy_factor);
        out<<rho<<','<<factor<<','<<std::log(Lambda)<<','<<c.Phi0_norm_bound<<','<<c.u0_norm_bound<<','
           <<c.ball_radius<<','<<c.Phi_domain_bound<<','<<c.u_domain_bound<<','
           <<c.source.phi_source<<','<<c.source.u_source<<','<<c.source.source_norm<<','
           <<c.lipschitz.contraction_bound<<','<<c.invariant_lhs<<','<<c.invariant_margin<<','
           <<(c.infinite_alpha_beta?1:0)<<','<<(c.inverse_certified?1:0)<<','
           <<(c.contraction_certified?1:0)<<','<<(c.invariant_ball_certified?1:0)<<','
           <<(c.proposition_b2_certificate?1:0)<<'\n';
    }
    return 0;
}
