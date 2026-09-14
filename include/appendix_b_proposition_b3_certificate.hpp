#pragma once

#include "appendix_b_axis_asymptotic.hpp"
#include "appendix_b_proposition_b2_certificate.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

// Point-evaluation constants implied directly by the B.4 alpha weights at beta=0:
// |F_alpha(eta)| <= ||F||_Brho 20^{-alpha}/(alpha+1)^2.
// Hence, for t=Y/20<1,
//   |F(Y,eta)| <= ||F|| S0(t),
//   |D_X F(Y,eta)| <= ||F|| S1(t), D_X=Y d_Y.
inline double appendix_b_brho_point_eval_S0(double Y) {
    if(!(Y>=0.0 && Y<20.0)) throw std::invalid_argument("B-rho point evaluation needs 0<=Y<20");
    const double t=Y/20.0;
    double s=0.0, p=1.0;
    for(std::size_t a=0;a<100000;++a){
        const double ap1=static_cast<double>(a+1);
        const double q=p/(ap1*ap1);
        s+=q;
        if(a>20 && std::abs(q)<1e-16*std::max(1.0,std::abs(s))) break;
        p*=t;
    }
    return s;
}

inline double appendix_b_brho_point_eval_S1(double Y) {
    if(!(Y>=0.0 && Y<20.0)) throw std::invalid_argument("B-rho point evaluation needs 0<=Y<20");
    const double t=Y/20.0;
    double s=0.0, p=1.0;
    for(std::size_t a=0;a<100000;++a){
        const double ap1=static_cast<double>(a+1);
        const double q=static_cast<double>(a)*p/(ap1*ap1);
        s+=q;
        if(a>20 && std::abs(q)<1e-16*std::max(1.0,std::abs(s))) break;
        p*=t;
    }
    return s;
}

struct AppendixBPropositionB3Certificate {
    double Ymax{};
    double Yexit{};
    double correction_brho_radius{};
    double point_eval_S0_Ymax{};
    double point_eval_S0_exit{};
    double point_eval_S1_exit{};
    double phi0_global_lower{};
    double phi_correction_sup{};
    double phi_lower{};
    double endpoint_phi_correction_sup{};
    double endpoint_DXphi_correction_sup{};
    double unperturbed_azimuthal_shear_lower{};
    double perturbed_azimuthal_shear_lower{};
    double endpoint_cex{};
    double log_C_chosen{};
    double log_C_required_axial{};
    double log_C_axial_margin{};
    bool inherited_from_b2{};
    bool phi_positive{};
    bool azimuthal_branch_certified{};
    bool axial_branch_certified{};
    bool endpoint_B19_certified{};
    bool proposition_b3_certificate{};
};

inline AppendixBPropositionB3Certificate appendix_b_proposition_b3_certificate(
    const AppendixBPropositionB2Certificate& b2,
    const AppendixBEndpointAudit& endpoint,
    double log_C_chosen,
    double Ymax=4.1,double Yexit=4.0,
    double phi0_lower=0.265,
    double target_exit=2.0) {
    if(!(Ymax>=Yexit && Yexit>0.0 && Ymax<20.0) || !(phi0_lower>0.0))
        throw std::invalid_argument("invalid Proposition B.3 continuation parameters");

    const double r=b2.ball_radius;
    const double S0max=appendix_b_brho_point_eval_S0(Ymax);
    const double S0=appendix_b_brho_point_eval_S0(Yexit);
    const double S1=appendix_b_brho_point_eval_S1(Yexit);
    const double e0max=r*S0max;
    const double e0=r*S0;
    const double e1=r*S1;
    const double philow=phi0_lower-e0max;

    // On chi>.99 the unperturbed endpoint shear is
    // p1=-2 D_X Phi0/Phi0 > 2.36.  The same coefficient norm controls both
    // Phi-Phi0 and D_X(Phi-Phi0).  Since the quotient lower bound is increasing
    // in the positive base Phi0, phi0_lower gives a uniform continuation bound.
    const double p0=endpoint.min_azimuthal_p1_on_chi_branch;
    double paz=-std::numeric_limits<double>::infinity();
    if(phi0_lower+e0>0.0)
        paz=(p0*phi0_lower-2.0*e1)/(phi0_lower+e0);

    // On the complementary branch the manuscript uses |Z*| separation and C
    // to make the axial term dominate.  The existing endpoint audit stores the
    // unperturbed logarithmic C threshold.  A relative Phi perturbation <=eps
    // is absorbed by increasing C by 1/(1-eps); this is conservative and tends
    // continuously to the asymptotic threshold as r->0.
    const double rel=(phi0_lower>0.0)?e0/phi0_lower:std::numeric_limits<double>::infinity();
    double log_required=std::numeric_limits<double>::infinity();
    if(rel<1.0)
        log_required=endpoint.max_log_required_C_normalized-std::log1p(-rel);
    const double logmargin=log_C_chosen-log_required;

    const bool inherited=b2.proposition_b2_certificate;
    const bool pos=inherited && philow>0.0;
    const bool az=inherited && pos && endpoint.azimuthal_branch_certified && paz>target_exit;
    const bool ax=inherited && pos && endpoint.axial_branch_separated && std::isfinite(logmargin) && logmargin>0.0;
    const bool b19=az && ax;
    const double cex=b19?std::min(paz-target_exit,0.3):0.0;
    const bool full=inherited && pos && b19;

    return {Ymax,Yexit,r,S0max,S0,S1,phi0_lower,e0max,philow,e0,e1,p0,paz,cex,
        log_C_chosen,log_required,logmargin,inherited,pos,az,ax,b19,full};
}

} // namespace nsblowup
