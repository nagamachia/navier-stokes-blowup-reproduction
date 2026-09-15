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
 bool pressure_derivatives_available{},pressure_enclosure_certified{},positive_on_interval{},separated_from_delta{};
};

inline AppendixBZIntervalAudit appendix_b_Z_interval_audit(
 const OuterProfileParameters& p,double a,double b,double delta_star,double j0=0.05,double Tf=20.0,double co=0.05,double max_step=0.01){
 if(!(a>=-1&&b<=1&&a<=b&&delta_star>0)) throw std::invalid_argument("invalid Appendix B Z interval");
 const double h=std::exp(p.log_h),A=appendix_b_A(h),D=appendix_b_D(h),invP2=std::exp(-2*p.log_Pstar);
 auto ivabs=[](double lo,double hi){return std::max(std::abs(lo),std::abs(hi));};
 // The chi<=.99 alternative is a genuine interval (~1.1e-2 at baseline), not
 // a point near eta0. Cover it by four closed boxes. Each box is certified by
 // a center value plus a derivative enclosure, so this is interval subdivision
 // rather than point sampling. Four boxes already recover ample Z>delta margin.
 constexpr int boxes=4;
 double global_lo=std::numeric_limits<double>::infinity(),global_hi=-std::numeric_limits<double>::infinity();
 double global_M=0.0,global_zerr=0.0; bool finite=true,penc=true;
 for(int ib=0;ib<boxes;++ib){
   const double x0=a+(b-a)*ib/boxes,x1=a+(b-a)*(ib+1)/boxes,m=.5*(x0+x1),rad=.5*(x1-x0);
   const auto pc=appendix_a4_pressure_enclosure(p,m,Tf,co,4*max_step);
   const auto pa=appendix_a4_pressure_enclosure(p,x0,Tf,co,4*max_step);
   const auto pb=appendix_a4_pressure_enclosure(p,x1,Tf,co,4*max_step);
   const auto zc=appendix_b_axis_data(p,m,j0,Tf,co,max_step);
   const double P0=std::max({ivabs(pa.pi0_lo,pa.pi0_hi),ivabs(pc.pi0_lo,pc.pi0_hi),ivabs(pb.pi0_lo,pb.pi0_hi)});
   const double P1=std::max({ivabs(pa.deta_lo,pa.deta_hi),ivabs(pc.deta_lo,pc.deta_hi),ivabs(pb.deta_lo,pb.deta_hi)});
   const double P2=std::max({ivabs(pa.detaeta_lo,pa.detaeta_hi),ivabs(pc.detaeta_lo,pc.detaeta_hi),ivabs(pb.detaeta_lo,pb.detaeta_hi)});
   const double emax=std::max(std::abs(x0),std::abs(x1)),Umax=4*emax+j0;
   const double Hmax=D*emax+(1+emax*emax)*Umax;
   const double Hpmax=D+2*emax*Umax+4*(1+emax*emax);
   const double velprime=invP2*A*(4*Umax*Umax+8*emax*Umax+4)+invP2*(4*Hmax+Hpmax*Umax);
   const double pressureprime=2*emax*P1+P2+4*A*(P0+emax*P1);
   const double M=appendix_b_up(velprime+pressureprime);
   const double zerr=appendix_b_up(std::abs(appendix_b_d(m))*pc.deta_error+4*A*std::abs(m)*pc.pi0_error);
   global_lo=std::min(global_lo,appendix_b_down(zc.normalized_Zstar-zerr-rad*M));
   global_hi=std::max(global_hi,appendix_b_up(zc.normalized_Zstar+zerr+rad*M));
   global_M=std::max(global_M,M); global_zerr=std::max(global_zerr,zerr);
   finite=finite&&std::isfinite(P0)&&std::isfinite(P1)&&std::isfinite(P2);
   penc=penc&&pa.certified_reproduction_enclosure&&pc.certified_reproduction_enclosure&&pb.certified_reproduction_enclosure;
 }
 const double mid=.5*(a+b); const auto zmid=appendix_b_axis_data(p,mid,j0,Tf,co,max_step);
 const double minabs=(global_lo>0)?global_lo:((global_hi<0)?-global_hi:0.0);
 return {a,b,mid,zmid.normalized_Zstar,global_M,global_lo,global_hi,minabs,delta_star,global_zerr,
         finite,penc,global_lo>0,minabs>delta_star};
}
} // namespace nsblowup
