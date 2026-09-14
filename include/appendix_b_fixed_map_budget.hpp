#pragma once

#include "appendix_b_brho.hpp"
#include "appendix_b_full_inverse.hpp"
#include "appendix_b_remainder_budget.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct AppendixBFixedMapBudget {
    double T_bound{}; // crude one-step bound, diagnostic only
    double inverse_one_plus_T_bound{}; // finite prefix + analytic factorial tail
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

    const auto na=op.radial_order, nb=op.eta_order;
    const double rho=op.rho, C=op.product;
    const double D=0.5-h;

    const double AXprod2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,false,false,true);
    const double AXdx2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,true,false,true);
    const double AXdeta2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,true,false,false,true);
    const double AXdx1=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,true,false,true);

    // The crude one-step norm can be >1.  The inverse instead uses an exact
    // finite Neumann prefix and the analytic factorial radial tail.
    const double t=0.5*op.product_J2*m.chi;
    const auto inv_audit=appendix_b_full_inverse_audit(na,nb,rho,m.chi,64,256);
    const double inv=inv_audit.full_inverse_bound;
    const bool inv_ok=inv_audit.tail_certified && std::isfinite(inv) && inv>0.0;

    const double outer1=inv_ok ? inv*C*m.invL/(2.0*Lambda) : INFINITY;
    const double outer2=C*m.invL/(2.0*Lambda);

    AppendixBFixedMapBudget out;
    out.T_bound=t;
    out.inverse_one_plus_T_bound=inv;
    out.inverse_finite_prefix_bound=inv_audit.finite_prefix_bound;
    out.inverse_analytic_tail_bound=inv_audit.analytic_tail_bound;
    out.inverse_factorial_K=inv_audit.factorial_envelope_K;
    out.inverse_envelope_ratio=inv_audit.max_envelope_ratio;
    out.inverse_certified=inv_ok;

    const double pair_uPhi=Phi_bound+u_bound;

    // ---- Phi map: (1+T)^-1 J2 R1 /(2 Lambda) ----
    out.phi_Wstar_Phi=outer1*op.product_J2*m.Wstar;
    out.phi_h_background_Phi=outer1*op.product_J2*h*m.one_minus_2etaUstar;

    // B=-2D eta A_X(u)-d d_eta A_X(u).
    out.phi_B_etaAX_u_Phi=outer1*(2.0*D/Lambda)*C*m.eta*AXprod2*pair_uPhi;
    out.phi_B_detaAX_u_Phi=outer1*(1.0/Lambda)*C*m.d*AXdeta2*pair_uPhi;
    out.phi_h_u_Phi=outer1*(2.0*h/Lambda)*C*m.eta*op.product_J2*pair_uPhi;
    out.phi_zeta_u_Phi=outer1*C*m.d*m.zeta*op.product_J2*pair_uPhi;

    out.phi_Wstar_DXPhi=outer1*op.DX_J2*m.Wstar;
    out.phi_B_etaAX_u_DXPhi=outer1*(2.0*D/Lambda)*C*m.eta*AXdx2*pair_uPhi;
    // Direct B.6: J2[(d_eta A_X u)(D_X Phi)].
    out.phi_B_detaAX_u_DXPhi=outer1*(1.0/Lambda)*C*m.d*op.mixed_AX_J2*pair_uPhi;

    out.phi_Hstar_detaPhi=outer1*op.deta_J2*m.Hstar;
    out.phi_u_detaPhi=outer1*(1.0/Lambda)*C*m.d*op.deta_J2*pair_uPhi;

    out.M_phi_map=out.phi_Wstar_Phi+out.phi_h_background_Phi+
        out.phi_B_etaAX_u_Phi+out.phi_B_detaAX_u_Phi+out.phi_h_u_Phi+
        out.phi_zeta_u_Phi+out.phi_Wstar_DXPhi+out.phi_B_etaAX_u_DXPhi+
        out.phi_B_detaAX_u_DXPhi+out.phi_Hstar_detaPhi+out.phi_u_detaPhi;

    // ---- u map: J1 R2 /(2 Lambda) ----
    out.u_linear=outer2*op.product_J1*m.linear_u;
    out.u_quadratic=outer2*(2.0*(0.5+h)/Lambda)*C*m.eta*op.product_J1*(2.0*u_bound);
    out.u_Wstar_DXu=outer2*op.DX_J1*m.Wstar;
    out.u_B_etaAX_u_DXu=outer2*(2.0*D/Lambda)*C*m.eta*AXdx1*(2.0*u_bound);
    // Direct B.6: J1[(d_eta A_X u)(D_X u)].
    out.u_B_detaAX_u_DXu=outer2*(1.0/Lambda)*C*m.d*op.mixed_AX_J1*(2.0*u_bound);
    out.u_Hstar_detau=outer2*op.deta_J1*m.Hstar;
    out.u_u_detau=outer2*(1.0/Lambda)*C*m.d*op.deta_J1*(2.0*u_bound);

    // Pressure terms remain conservative; they are not the mixed B.6 products.
    const double Lp=2.0*op.I*C*C*C*g_bound*g_bound*Phi_bound;
    out.u_pressure_p=outer2*op.J1*C*m.Aeta4*Lp;
    out.u_pressure_peta=outer2*op.J1*C*m.d*op.deta*Lp;
    out.u_pressure_DXp=outer2*op.J1*C*(2.0*m.eta)*op.DX*Lp;

    out.M_u_map=out.u_linear+out.u_quadratic+out.u_Wstar_DXu+
        out.u_B_etaAX_u_DXu+out.u_B_detaAX_u_DXu+out.u_Hstar_detau+
        out.u_u_detau+out.u_pressure_p+out.u_pressure_peta+out.u_pressure_DXp;
    out.contraction_bound=std::max(out.M_phi_map,out.M_u_map);
    out.contraction_certified=inv_ok && std::isfinite(out.contraction_bound) && out.contraction_bound<1.0;
    return out;
}

} // namespace nsblowup
