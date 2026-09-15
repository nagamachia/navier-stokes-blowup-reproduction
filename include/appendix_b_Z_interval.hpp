#pragma once

#include "appendix_b_endpoint_interval.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBZIntervalAudit {
    double eta_lo{};
    double eta_hi{};
    double eta_center{};
    double Z_center{};
    double max_abs_Zprime{};
    double Z_lower{};
    double Z_upper{};
    double min_abs_Z{};
    double delta_star{};
    bool pressure_derivatives_available{};
    bool positive_on_interval{};
    bool separated_from_delta{};
};

// Rigorous-in-eta enclosure for the narrow axial candidate interval.
// Pressure values and first two eta derivatives are evaluated at the center;
// the derivative envelope is built from the explicit B.1 formula and endpoint
// samples of the smooth A.4 datum, enlarged by nextafter rounding.  This removes
// the old uniform eta grid from the B.19 axial decision.  The pressure evaluator
// itself is still floating-point RK/Simpson, so this is a reproduction interval
// certificate rather than a formal interval integration of Appendix A.4.
inline AppendixBZIntervalAudit appendix_b_Z_interval_audit(
    const OuterProfileParameters& p,double a,double b,double delta_star,
    double j0=0.05,double Tf=20.0,double co=0.05,double max_step=0.01){
    if(!(a>=-1.0&&b<=1.0&&a<=b&&delta_star>0.0))
        throw std::invalid_argument("invalid Appendix B Z interval");
    const double m=0.5*(a+b), rad=0.5*(b-a);
    const double h=std::exp(p.log_h), A=appendix_b_A(h), D=appendix_b_D(h);
    const double invP2=std::exp(-2.0*p.log_Pstar);
    const auto pc=appendix_a4_pressure_datum(p,m,Tf,co,max_step);
    const auto pa=appendix_a4_pressure_datum(p,a,Tf,co,max_step);
    const auto pb=appendix_a4_pressure_datum(p,b,Tf,co,max_step);
    const auto zc=appendix_b_axis_data(p,m,j0,Tf,co,max_step);

    auto absmax=[](double x,double y,double z){return std::max({std::abs(x),std::abs(y),std::abs(z)});};
    const double P0=absmax(pa.normalized_pi0,pc.normalized_pi0,pb.normalized_pi0);
    const double P1=absmax(pa.normalized_deta,pc.normalized_deta,pb.normalized_deta);
    const double P2=absmax(pa.normalized_detaeta,pc.normalized_detaeta,pb.normalized_detaeta);

    // Differentiate B.1.  U=4 eta+j0, H=D eta+(1-eta^2)U.
    // The velocity contribution is negligible after P_*^{-2}, but is retained.
    const double emax=std::max(std::abs(a),std::abs(b));
    const double Umax=4.0*emax+j0;
    const double Hmax=D*emax+(1.0+emax*emax)*Umax;
    const double Hpmax=D+2.0*emax*Umax+4.0*(1.0+emax*emax);
    const double velprime=invP2*A*(4.0*Umax*Umax + 8.0*emax*Umax + 4.0)
                         + invP2*(4.0*Hmax + Hpmax*Umax);
    // d/deta[-d Pi' + 4 A eta Pi] = 2 eta Pi' - d Pi'' + 4 A(Pi+eta Pi').
    const double pressureprime=2.0*emax*P1 + P2 + 4.0*A*(P0+emax*P1);
    const double M=appendix_b_up(velprime+pressureprime);
    const double lo=appendix_b_down(zc.normalized_Zstar-rad*M);
    const double hi=appendix_b_up(zc.normalized_Zstar+rad*M);
    const double minabs=(lo>0.0)?lo:((hi<0.0)?-hi:0.0);
    return {a,b,m,zc.normalized_Zstar,M,lo,hi,minabs,delta_star,
            std::isfinite(P0)&&std::isfinite(P1)&&std::isfinite(P2),
            lo>0.0,minabs>delta_star};
}

} // namespace nsblowup
