#pragma once

#include "appendix_b_brho.hpp"
#include "appendix_b_full_inverse.hpp"
#include "appendix_b_remainder_budget.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace nsblowup {

struct AppendixBFixedMapBudget {
    double T_bound{};
    double inverse_one_plus_T_bound{};
    std::size_t inverse_finite_k{};
    double inverse_finite_prefix_bound{};
    double inverse_analytic_tail_bound{};
    double inverse_factorial_K{};
    double inverse_envelope_ratio{};

    double phi_Wstar_Phi{};
    double phi_h_background_Phi{};
    double phi_B_etaAX_u_Phi{};
    double phi_B_detaAX_u_Phi{};
    double phi_h_u_Phi{};
    double phi_zeta_u_Phi{};
    double phi_Wstar_DXPhi{};
    double phi_B_etaAX_u_DXPhi{};
    double phi_B_detaAX_u_DXPhi{};
    double phi_Hstar_detaPhi{};
    double phi_u_detaPhi{};

    double u_linear{};
    double u_quadratic{};
    double u_Wstar_DXu{};
    double u_B_etaAX_u_DXu{};
    double u_B_detaAX_u_DXu{};
    double u_Hstar_detau{};
    double u_u_detau{};
    double u_pressure_p{};
    double u_pressure_peta{};
    double u_pressure_DXp{};

    double M_phi_map{};
    double M_u_map{};
    double contraction_bound{};
    bool inverse_certified{};
    bool infinite_eta_certified{};
    bool infinite_radial_certified{};
    bool contraction_certified{};
};

inline AppendixBFixedMapBudget appendix_b_fixed_map_lipschitz_budget(
    const AppendixBRhoOperatorAudit& op,
    const AppendixBMultiplierNorms& m,
    double h,
    double Lambda,
    double Phi_bound,
    double u_bound,
    double g_bound) {
    if (!(h>0.0 && h<0.5) || !(Lambda>0.0) || !(Phi_bound>0.0) ||
        !(u_bound>=0.0) || !(g_bound>=0.0) || !(m.invL>0.0) || !(m.chi>0.0))
        throw std::invalid_argument("invalid Appendix B fixed-map budget parameters");

    const double C=op.product;
    const double D=0.5-h;
    const double AXprod2=op.AX_product_J2;
    const double AXdx2=op.AX_DX_J2;
    const double AXdeta2=op.AX_deta_J2;
    const double AXdx1=op.AX_DX_J1;

    const double t=0.5*op.product_J2*m.chi;
    AppendixBFullInverseAudit inv_audit;
    if(op.infinite_eta_certified && op.infinite_radial_certified) {
        // The factorial envelope already controls every starting radial degree,
        // so the prefix length is only a summation split, not an alpha cutoff.
        inv_audit=appendix_b_full_inverse_all_beta_audit(op.rho,m.chi,1);
    } else if(op.infinite_eta_certified) {
        inv_audit=appendix_b_full_inverse_all_beta_audit(op.rho,m.chi,op.radial_order);
    } else {
        const double K=appendix_b_T_factorial_envelope_K(m.chi);
        std::size_t finite_k=op.radial_order;
        while(K/(static_cast<double>(finite_k+2)*static_cast<double>(finite_k+3))>=0.5){
            ++finite_k;
            if(finite_k>160)
                throw std::runtime_error("Appendix B inverse prefix exceeds safe direct-double range");
        }
        inv_audit=appendix_b_full_inverse_audit(
            finite_k,op.eta_order,op.rho,m.chi,64,128);
    }
    const double inv=inv_audit.full_inverse_bound;
    const bool inv_ok=inv_audit.tail_certified && std::isfinite(inv) && inv>0.0;

    const double outer1=inv_ok ? inv*C*m.invL/(2.0*Lambda) : INFINITY;
    const double outer2=C*m.invL/(2.0*Lambda);

    AppendixBFixedMapBudget out;
    out.T_bound=t;
    out.inverse_one_plus_T_bound=inv;
    out.inverse_finite_k=inv_audit.finite_k;
    out.inverse_finite_prefix_bound=inv_audit.finite_prefix_bound;
    out.inverse_analytic_tail_bound=inv_audit.analytic_tail_bound;
    out.inverse_factorial_K=inv_audit.factorial_envelope_K;
    out.inverse_envelope_ratio=inv_audit.max_envelope_ratio;
    out.inverse_certified=inv_ok;
    out.infinite_eta_certified=op.infinite_eta_certified;
    out.infinite_radial_certified=op.infinite_radial_certified;

    const double pair_uPhi=Phi_bound+u_bound;

    out.phi_Wstar_Phi=outer1*op.product_J2*m.Wstar;
    out.phi_h_background_Phi=outer1*op.product_J2*h*m.one_minus_2etaUstar;
    out.phi_B_etaAX_u_Phi=outer1*(2.0*D/Lambda)*C*m.eta*AXprod2*pair_uPhi;
    out.phi_B_detaAX_u_Phi=outer1*(1.0/Lambda)*C*m.d*AXdeta2*pair_uPhi;
    out.phi_h_u_Phi=outer1*(2.0*h/Lambda)*C*m.eta*op.product_J2*pair_uPhi;
    out.phi_zeta_u_Phi=outer1*C*m.d*m.zeta*op.product_J2*pair_uPhi;
    out.phi_Wstar_DXPhi=outer1*op.DX_J2*m.Wstar;
    out.phi_B_etaAX_u_DXPhi=outer1*(2.0*D/Lambda)*C*m.eta*AXdx2*pair_uPhi;
    out.phi_B_detaAX_u_DXPhi=outer1*(1.0/Lambda)*C*m.d*op.mixed_AX_J2*pair_uPhi;
    out.phi_Hstar_detaPhi=outer1*op.deta_J2*m.Hstar;
    out.phi_u_detaPhi=outer1*(1.0/Lambda)*C*m.d*op.deta_J2*pair_uPhi;

    out.M_phi_map=out.phi_Wstar_Phi+out.phi_h_background_Phi+
        out.phi_B_etaAX_u_Phi+out.phi_B_detaAX_u_Phi+out.phi_h_u_Phi+
        out.phi_zeta_u_Phi+out.phi_Wstar_DXPhi+out.phi_B_etaAX_u_DXPhi+
        out.phi_B_detaAX_u_DXPhi+out.phi_Hstar_detaPhi+out.phi_u_detaPhi;

    out.u_linear=outer2*op.product_J1*m.linear_u;
    out.u_quadratic=outer2*(2.0*(0.5+h)/Lambda)*C*m.eta*op.product_J1*(2.0*u_bound);
    out.u_Wstar_DXu=outer2*op.DX_J1*m.Wstar;
    out.u_B_etaAX_u_DXu=outer2*(2.0*D/Lambda)*C*m.eta*AXdx1*(2.0*u_bound);
    out.u_B_detaAX_u_DXu=outer2*(1.0/Lambda)*C*m.d*op.mixed_AX_J1*(2.0*u_bound);
    out.u_Hstar_detau=outer2*op.deta_J1*m.Hstar;
    out.u_u_detau=outer2*(1.0/Lambda)*C*m.d*op.deta_J1*(2.0*u_bound);

    const double pressure_source_lip=2.0*C*C*C*g_bound*g_bound*Phi_bound;
    out.u_pressure_p=outer2*C*m.Aeta4*op.pressure_J1_I*pressure_source_lip;
    out.u_pressure_peta=outer2*C*m.d*op.pressure_J1_detaI*pressure_source_lip;
    out.u_pressure_DXp=outer2*C*(2.0*m.eta)*op.pressure_J1_DXI*pressure_source_lip;

    out.M_u_map=out.u_linear+out.u_quadratic+out.u_Wstar_DXu+
        out.u_B_etaAX_u_DXu+out.u_B_detaAX_u_DXu+out.u_Hstar_detau+
        out.u_u_detau+out.u_pressure_p+out.u_pressure_peta+out.u_pressure_DXp;
    out.contraction_bound=std::max(out.M_phi_map,out.M_u_map);
    out.contraction_certified=inv_ok && op.infinite_eta_certified &&
        op.infinite_radial_certified && std::isfinite(out.contraction_bound) &&
        out.contraction_bound<1.0;
    return out;
}

} // namespace nsblowup
