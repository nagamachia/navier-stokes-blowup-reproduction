#pragma once

#include "appendix_b_endpoint_interval.hpp"
#include "appendix_a4_pressure_enclosure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBZIntervalAudit {
    double eta_lo{},eta_hi{},eta_center{};
    double Z_center{},max_abs_Zprime{},Z_lower{},Z_upper{},min_abs_Z{},delta_star{};
    double pressure_center_error{};
    bool pressure_derivatives_available{};
    bool pressure_enclosure_certified{};
    bool positive_on_interval{};
    bool separated_from_delta{};
};

inline AppendixBZIntervalAudit appendix_b_Z_interval_audit(
    const OuterProfileParameters& p,double a,double b,double delta_star,
    double j0=0.05,double Tf=20.0,double co=0.05,double max_step=0.01){
    if(!(a>=-1.0&&b<=1.0&&a<=b&&delta_star>0.0)) throw std::invalid_argument("invalid Appendix B Z interval");
    const double m=0.5*(a+b),rad=0.5*(b-a);
    const double h=std::exp(p.log_h),A=appendix_b_A(h),D=appendix_b_D(h),invP2=std::exp(-2.0*p.log_Pstar);
    // Use three nested radial resolutions; coarse_step=4*max_step makes the
    // finest member exactly the caller's requested max_step.
    const auto pc=appendix_a4_pressure_enclosure(p,m,Tf,co,4.0*max_step);
    const auto pa=appendix_a4_pressure_enclosure(p,a,Tf,co,4.0*max_step);
    const auto pb=appendix_a4_pressure_enclosure(p,b,Tf,co,4.0*max_step);
    const auto zc=appendix_b_axis_data(p,m,j0,Tf,co,max_step);
    auto ivabs=[](double lo,double hi){return std::max(std::abs(lo),std::abs(hi));};
    const double P0=std::max({ivabs(pa.pi0_lo,pa.pi0_hi),ivabs(pc.pi0_lo,pc.pi0_hi),ivabs(pb.pi0_lo,pb.pi0_hi)});
    const double P1=std::max({ivabs(pa.deta_lo,pa.deta_hi),ivabs(pc.deta_lo,pc.deta_hi),ivabs(pb.deta_lo,pb.deta_hi)});
    const double P2=std::max({ivabs(pa.detaeta_lo,pa.detaeta_hi),ivabs(pc.detaeta_lo,pc.detaeta_hi),ivabs(pb.detaeta_lo,pb.detaeta_hi)});
    const double emax=std::max(std::abs(a),std::abs(b));
    const double Umax=4.0*emax+j0;
    const double Hmax=D*emax+(1.0+emax*emax)*Umax;
    const double Hpmax=D+2.0*emax*Umax+4.0*(1.0+emax*emax);
    const double velprime=invP2*A*(4.0*Umax*Umax+8.0*emax*Umax+4.0)+invP2*(4.0*Hmax+Hpmax*Umax);
    const double pressureprime=2.0*emax*P1+P2+4.0*A*(P0+emax*P1);
    const double M=appendix_b_up(velprime+pressureprime);
    // Center Z pressure uncertainty: |-d delta Pi' +4A eta delta Pi|.
    const double d=appendix_b_d(m);
    const double zerr=appendix_b_up(std::abs(d)*pc.deta_error+4.0*A*std::abs(m)*pc.pi0_error);
    const double lo=appendix_b_down(zc.normalized_Zstar-zerr-rad*M);
    const double hi=appendix_b_up(zc.normalized_Zstar+zerr+rad*M);
    const double minabs=(lo>0.0)?lo:((hi<0.0)?-hi:0.0);
    const bool penc=pa.certified_reproduction_enclosure&&pc.certified_reproduction_enclosure&&pb.certified_reproduction_enclosure;
    return {a,b,m,zc.normalized_Zstar,M,lo,hi,minabs,delta_star,zerr,
      std::isfinite(P0)&&std::isfinite(P1)&&std::isfinite(P2),penc,lo>0.0,minabs>delta_star};
}

} // namespace nsblowup
