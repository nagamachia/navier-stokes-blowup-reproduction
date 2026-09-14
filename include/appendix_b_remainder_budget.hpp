#pragma once

#include "appendix_b_brho.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct AppendixBMultiplierNorms {
    // B-rho norms of scalar eta-dependent multipliers in the exact B.15 remainder.
    double invL{};          // 1/L
    double Wstar{};
    double Hstar{};
    double eta{};
    double d{};
    double zeta{};
    double chi{};           // H_*^2/(H_*^2+sigma_*^2), for T=(1/2)J_2(chi .)
    double one_minus_2etaUstar{};
    double linear_u{};      // A(1-4 eta U*)+d U*_eta
    double Aeta4{};         // 4 A eta
};

struct AppendixBRemainderBudget {
    double r1_bracket_phi{};
    double r1_W_DXPhi{};
    double r1_H_detaPhi{};
    double r2_linear_u{};
    double r2_u2{};
    double r2_W_DXu{};
    double r2_H_deta_u{};
    double r2_u_deta_u{};
    double r2_pressure_p{};
    double r2_pressure_peta{};
    double r2_pressure_DXp{};
    double M_R1{};
    double M_R2{};
    double M{};
    double pressure_fraction{};
};

inline AppendixBRemainderBudget appendix_b_remainder_lipschitz_budget(
    const AppendixBRhoOperatorAudit& op,
    const AppendixBMultiplierNorms& m,
    double h,
    double Lambda,
    double Phi_bound,
    double u_bound,
    double g_bound) {
    if (!(h>0.0) || !(Lambda>0.0) || !(Phi_bound>0.0) ||
        !(u_bound>=0.0) || !(g_bound>=0.0) || !(m.invL>0.0))
        throw std::invalid_argument("invalid Appendix B remainder budget parameters");

    const double C=op.product, DX=op.DX, AX=op.AX, De=op.deta, I=op.I;
    const double LB = C*m.eta*AX + C*m.d*AX*De;
    const double Wb = m.Wstar + LB*u_bound/Lambda;
    const double LW = LB/Lambda;
    const double LH = C*m.d/Lambda;
    const double Hb = m.Hstar + LH*u_bound;

    const double Ku = C*(2.0*h*m.eta/Lambda + C*m.d*m.zeta);
    const double Kb = Wb + h*m.one_minus_2etaUstar + Ku*u_bound;
    const double LK = LW + Ku;

    AppendixBRemainderBudget out;
    out.r1_bracket_phi = C*m.invL * C*(Kb + LK*Phi_bound);
    out.r1_W_DXPhi = C*m.invL * C*DX*(Wb + LW*Phi_bound);
    out.r1_H_detaPhi = C*m.invL * C*De*(Hb + LH*Phi_bound);
    out.M_R1 = out.r1_bracket_phi + out.r1_W_DXPhi + out.r1_H_detaPhi;

    out.r2_linear_u = C*m.invL * C*m.linear_u;
    out.r2_u2 = C*m.invL * (C*C)*(2.0*(0.5+h)*m.eta/Lambda)*(2.0*u_bound);
    out.r2_W_DXu = C*m.invL * C*DX*(Wb + LW*u_bound);
    out.r2_H_deta_u = C*m.invL * C*m.Hstar*De;
    out.r2_u_deta_u = C*m.invL * (C*C)*(m.d/Lambda)*(2.0*De*u_bound);

    const double Lp = 2.0*I*C*C*C*g_bound*g_bound*Phi_bound;
    out.r2_pressure_p = C*m.invL * C*m.Aeta4*Lp;
    out.r2_pressure_peta = C*m.invL * C*m.d*De*Lp;
    out.r2_pressure_DXp = C*m.invL * C*(2.0*m.eta)*DX*Lp;

    out.M_R2 = out.r2_linear_u + out.r2_u2 + out.r2_W_DXu +
        out.r2_H_deta_u + out.r2_u_deta_u + out.r2_pressure_p +
        out.r2_pressure_peta + out.r2_pressure_DXp;
    out.M = std::max(out.M_R1,out.M_R2);
    const double pressure = out.r2_pressure_p + out.r2_pressure_peta + out.r2_pressure_DXp;
    out.pressure_fraction = out.M_R2>0.0 ? pressure/out.M_R2 : 0.0;
    return out;
}

} // namespace nsblowup
