#pragma once

#include "appendix_b_fixed_map_budget.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBFixedMapSourceBudget {
    double phi_source{};
    double u_source{};
    double source_norm{};
};

inline AppendixBFixedMapSourceBudget appendix_b_fixed_map_source_budget(
    const AppendixBRhoOperatorAudit& op,
    const AppendixBMultiplierNorms& m,
    double h,double Lambda,double Phi0,double u0,double g_bound) {
    if(!(h>0.0 && h<0.5) || !(Lambda>0.0) || !(Phi0>0.0) || !(u0>=0.0))
        throw std::invalid_argument("invalid Appendix B fixed-map source parameters");
    const double C=op.product;
    const double D=0.5-h;
    const double inv=appendix_b_full_inverse_all_beta_audit(op.rho,m.chi).full_inverse_bound;
    const double outer1=inv*C*m.invL/(2.0*Lambda);
    const double outer2=C*m.invL/(2.0*Lambda);

    double s1=0.0;
    s1 += outer1*op.product_J2*m.Wstar*Phi0;
    s1 += outer1*op.product_J2*h*m.one_minus_2etaUstar*Phi0;
    s1 += outer1*(2.0*D/Lambda)*C*m.eta*op.AX_product_J2*u0*Phi0;
    s1 += outer1*(1.0/Lambda)*C*m.d*op.AX_deta_J2*u0*Phi0;
    s1 += outer1*(2.0*h/Lambda)*C*m.eta*op.product_J2*u0*Phi0;
    s1 += outer1*C*m.d*m.zeta*op.product_J2*u0*Phi0;
    s1 += outer1*op.DX_J2*m.Wstar*Phi0;
    s1 += outer1*(2.0*D/Lambda)*C*m.eta*op.AX_DX_J2*u0*Phi0;
    s1 += outer1*(1.0/Lambda)*C*m.d*op.mixed_AX_J2*u0*Phi0;
    s1 += outer1*op.deta_J2*m.Hstar*Phi0;
    s1 += outer1*(1.0/Lambda)*C*m.d*op.deta_J2*u0*Phi0;

    double s2=0.0;
    s2 += outer2*op.product_J1*m.linear_u*u0;
    s2 += outer2*(2.0*(0.5+h)/Lambda)*C*m.eta*op.product_J1*u0*u0;
    s2 += outer2*op.DX_J1*m.Wstar*u0;
    s2 += outer2*(2.0*D/Lambda)*C*m.eta*op.AX_DX_J1*u0*u0;
    s2 += outer2*(1.0/Lambda)*C*m.d*op.mixed_AX_J1*u0*u0;
    s2 += outer2*op.deta_J1*m.Hstar*u0;
    s2 += outer2*(1.0/Lambda)*C*m.d*op.deta_J1*u0*u0;
    const double pressure_source=C*C*C*g_bound*g_bound*Phi0*Phi0;
    s2 += outer2*C*m.Aeta4*op.pressure_J1_I*pressure_source;
    s2 += outer2*C*m.d*op.pressure_J1_detaI*pressure_source;
    s2 += outer2*C*(2.0*m.eta)*op.pressure_J1_DXI*pressure_source;
    return {s1,s2,std::max(s1,s2)};
}

struct AppendixBPropositionB2Certificate {
    double Lambda{};
    double Phi0_norm_bound{};
    double u0_norm_bound{};
    double ball_radius{};
    double Phi_domain_bound{};
    double u_domain_bound{};
    AppendixBFixedMapSourceBudget source{};
    AppendixBFixedMapBudget lipschitz{};
    double invariant_lhs{};
    double invariant_margin{};
    bool infinite_alpha_beta{};
    bool inverse_certified{};
    bool contraction_certified{};
    bool invariant_ball_certified{};
    bool proposition_b2_certificate{};
};

inline AppendixBPropositionB2Certificate appendix_b_proposition_b2_certificate(
    const AppendixBRhoOperatorAudit& op,const AppendixBMultiplierNorms& m,
    double h,double Lambda,double u0_norm_bound,double g_bound=1.0) {
    if(!op.infinite_eta_certified || !op.infinite_alpha_certified)
        throw std::invalid_argument("Proposition B.2 certificate requires infinite-alpha-beta operators");
    const auto ia=appendix_b_full_inverse_all_beta_audit(op.rho,m.chi);
    const double Phi0=ia.full_inverse_bound; // ||(1+T)^-1 1|| <= inverse norm.
    double r=1.0;
    AppendixBFixedMapSourceBudget src{};
    AppendixBFixedMapBudget lip{};
    for(int it=0;it<12;++it){
        const double Pb=Phi0+r;
        const double ub=u0_norm_bound+r;
        lip=appendix_b_fixed_map_lipschitz_budget(op,m,h,Lambda,Pb,ub,g_bound);
        src=appendix_b_fixed_map_source_budget(op,m,h,Lambda,Phi0,u0_norm_bound,g_bound);
        if(!(lip.contraction_bound<1.0) || !std::isfinite(src.source_norm)) break;
        const double need=src.source_norm/(1.0-lip.contraction_bound);
        const double next=1.05*std::max(need,1e-300);
        if(next<=r) { r=std::max(next,1e-300); break; }
        r=next;
    }
    const double Pb=Phi0+r;
    const double ub=u0_norm_bound+r;
    lip=appendix_b_fixed_map_lipschitz_budget(op,m,h,Lambda,Pb,ub,g_bound);
    src=appendix_b_fixed_map_source_budget(op,m,h,Lambda,Phi0,u0_norm_bound,g_bound);
    const double lhs=src.source_norm+lip.contraction_bound*r;
    const double margin=r-lhs;
    const bool invball=std::isfinite(lhs) && lhs<=r && lip.contraction_bound<1.0;
    const bool full=op.infinite_eta_certified && op.infinite_alpha_certified &&
        lip.inverse_certified && lip.contraction_certified && invball;
    return {Lambda,Phi0,u0_norm_bound,r,Pb,ub,src,lip,lhs,margin,
        true,lip.inverse_certified,lip.contraction_certified,invball,full};
}

} // namespace nsblowup
