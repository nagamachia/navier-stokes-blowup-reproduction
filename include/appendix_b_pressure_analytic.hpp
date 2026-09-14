#pragma once

#include "appendix_a4_pressure.hpp"
#include "appendix_b_analytic_multiplier.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBPressureAnalyticAudit {
    double rho{};
    double inner_radius{};
    double outer_radius{};
    double f_sup_outer{};
    double pi0_normalized_sup_outer{};
    double dpi0_normalized_sup_inner{};
    double pi0_brho_norm{};
    double dpi0_brho_norm{};
    double cauchy_factor_pi{};
    double cauchy_factor_dpi{};
    bool pole_free{};
    bool certified{};
};

inline double appendix_b_cauchy_first_derivative_brho_factor(double rho,double R) {
    if(!(rho>0.0 && R>rho)) throw std::invalid_argument("first-derivative Cauchy radius must exceed rho");
    const double q=rho/R;
    double best=1.0/R;
    for(std::size_t b=0;b<100000;++b){
        // |d^(b+1)f| <= M (b+1)! R^(-b-1). Divide by the B.4 alpha=0
        // weight rho^(-b) b!/(b+1)^2.
        const double t=(1.0/R)*std::pow(q,static_cast<double>(b))*
            std::pow(static_cast<double>(b+1),3.0);
        best=std::max(best,t);
        const double ratio=q*std::pow(static_cast<double>(b+2)/static_cast<double>(b+1),3.0);
        if(b>12 && ratio<1.0 && t<best*1e-15) break;
    }
    return best;
}

// Analytic majorant for Appendix A.4 pressure datum.  Before A.10 every eta
// dependent contribution is K_pre f(eta)^2. During A.10 the integrand has the
// form c(y) f(eta)^(2 theta(y)), 0<=theta<=1.  At the A.10 endpoint the eta
// dependence cancels exactly, hence all later exterior/terminal contributions
// are constants in eta.  We recover the positive constants K_pre and K_post
// from the real-axis datum at eta=0, where f=1, and majorize the A.10 integral
// by its eta=0 value times max(1,|f|^2) on the complex stadium.
inline AppendixBPressureAnalyticAudit appendix_b_pressure_analytic_audit(
    const OuterProfileParameters& p,double rho,double inner_radius,double outer_radius,
    double Tf=20.0,double co=0.05,double max_step=0.02) {
    if(!(rho>0.0 && inner_radius>rho && outer_radius>inner_radius && outer_radius<1.0))
        throw std::invalid_argument("invalid Appendix B pressure analytic radii");

    // f(z)=1/[(z-i)(z+i)]. For dist(z,[-1,1])<=R<1, both factors are at
    // least 1-R from their poles, so |f| <= (1-R)^-2.
    const double f_sup=1.0/std::pow(1.0-outer_radius,2.0);

    // Reconstruct the eta-dependent positive radial masses from the same A.4
    // schedule at eta=0.  At eta=0, d_eta terms vanish and f=1.
    const auto d0=appendix_a4_pressure_datum(p,0.0,Tf,co,max_step);
    const double total_mass0=-2.0*d0.normalized_pi0;
    if(!(total_mass0>0.0)) throw std::runtime_error("Appendix A.4 pressure mass must be positive");

    // The inner/pre-A10 mass is at least the exact inner contribution 5.  For a
    // rigorous reproduction majorant we may safely assign the whole positive
    // eta=0 mass the worst f^2 amplification. This deliberately overestimates
    // A.10 and post-A10 constants but removes any need to complexify the radial
    // RK schedule.
    const double pi_sup=0.5*total_mass0*std::max(1.0,f_sup*f_sup);

    // First derivative on the inner stadium from Cauchy on the outer stadium.
    const double gap=outer_radius-inner_radius;
    const double dpi_sup=pi_sup/gap;

    // Full B-rho norms on the real eta segment. Pi uses the outer radius; Pi_eta
    // uses the available gap after first differentiating.
    const double Fpi=appendix_b_cauchy_brho_factor(rho,outer_radius);
    const double Fdpi=appendix_b_cauchy_first_derivative_brho_factor(rho,outer_radius);
    const double pi_norm=pi_sup*Fpi;
    const double dpi_norm=pi_sup*Fdpi;

    AppendixBPressureAnalyticAudit out;
    out.rho=rho;
    out.inner_radius=inner_radius;
    out.outer_radius=outer_radius;
    out.f_sup_outer=f_sup;
    out.pi0_normalized_sup_outer=pi_sup;
    out.dpi0_normalized_sup_inner=dpi_sup;
    out.pi0_brho_norm=pi_norm;
    out.dpi0_brho_norm=dpi_norm;
    out.cauchy_factor_pi=Fpi;
    out.cauchy_factor_dpi=Fdpi;
    out.pole_free=outer_radius<1.0;
    out.certified=out.pole_free && std::isfinite(pi_norm) && std::isfinite(dpi_norm) &&
        pi_norm>0.0 && dpi_norm>0.0;
    return out;
}

struct AppendixBZoverLAnalyticAudit {
    AppendixBPressureAnalyticAudit pressure;
    double normalized_Z_over_L_brho{}; // bound divided by P_*^2
    double log_Z_over_L_brho{};        // physical Z*/L bound in log scale
    double log_u0_brho_Ymax{};
    bool certified{};
};

inline AppendixBZoverLAnalyticAudit appendix_b_Z_over_L_analytic_audit(
    const OuterProfileParameters& p,double j0,double rho,double inner_radius,double outer_radius,
    double Ymax=4.1,double Tf=20.0,double co=0.05,double max_step=0.02) {
    const auto pa=appendix_b_pressure_analytic_audit(p,rho,inner_radius,outer_radius,Tf,co,max_step);
    if(!pa.certified) throw std::runtime_error("pressure analytic audit failed");
    const double h=std::exp(p.log_h);
    const double A=appendix_b_A(h);
    const double F=appendix_b_cauchy_brho_factor(rho,inner_radius);
    const double zmax=1.0+inner_radius;
    const double invLsup=1.0/(1.0-2.0*h*zmax*zmax);
    const double Usup=4.0*zmax+j0;
    const double Hsup=appendix_b_poly_sup_from_coeffs({j0,4.5-h,-j0,-4.0},inner_radius);
    const double dsup=1.0+zmax*zmax;

    // Z*/P_*^2 = P_*^-2[-A(1-2 eta U*)U* - H* U*_eta]
    //             - d (Pi0_eta/P_*^2) + 4 A eta (Pi0/P_*^2).
    const double invP2=std::exp(-2.0*p.log_Pstar);
    const double velocity_sup=invP2*(A*(1.0+2.0*zmax*Usup)*Usup+4.0*Hsup);
    const double pressure_sup=dsup*pa.dpi0_normalized_sup_inner+
        4.0*A*zmax*pa.pi0_normalized_sup_outer;
    const double znorm_sup=invLsup*(velocity_sup+pressure_sup);
    const double znorm_brho=F*znorm_sup;
    const double logZ=2.0*p.log_Pstar+std::log(znorm_brho);
    const double logu=std::log(Ymax/2.0)+logZ;
    return {pa,znorm_brho,logZ,logu,pa.certified && std::isfinite(logu)};
}

} // namespace nsblowup
