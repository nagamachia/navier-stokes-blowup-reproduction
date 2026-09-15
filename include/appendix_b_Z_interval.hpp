#pragma once
#include "appendix_b_endpoint_interval.hpp"
#include "appendix_a4_pressure_enclosure.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace nsblowup {
struct AppendixBZIntervalAudit { double eta_lo{},eta_hi{},eta_center{},Z_center{},max_abs_Zprime{},Z_lower{},Z_upper{},min_abs_Z{},delta_star{},pressure_center_error{}; bool pressure_derivatives_available{},pressure_enclosure_certified{},positive_on_interval{},separated_from_delta{}; };
inline AppendixBZIntervalAudit appendix_b_Z_interval_audit(const OuterProfileParameters& p,double a,double b,double delta_star,double j0=0.05,double Tf=20.0,double co=0.05,double max_step=0.01){
 if(!(a>=-1&&b<=1&&a<=b&&delta_star>0))throw std::invalid_argument("invalid Appendix B Z interval");
 const double h=std::exp(p.log_h),A=appendix_b_A(h),D=appendix_b_D(h),invP2=std::exp(-2*p.log_Pstar);
 auto ivabs=[](double lo,double hi){return std::max(std::abs(lo),std::abs(hi));};
 constexpr int boxes=4; double glo=std::numeric_limits<double>::infinity(),ghi=-glo,gM=0,gerr=0; bool finite=true,penc=true;
 for(int ib=0;ib<boxes;++ib){
  const double x0=a+(b-a)*ib/boxes,x1=a+(b-a)*(ib+1)/boxes,m=.5*(x0+x1),rad=.5*(x1-x0);
  // enclosure internally evaluates at coarse_step/4.  Passing max_step here
  // therefore makes the analytic h^4 truncation radius use max_step/4, while
  // the independent center value below remains at max_step.  This is not
  // Richardson: no resolution difference enters the error estimate.
  const auto pc=appendix_a4_pressure_enclosure(p,m,Tf,co,max_step),pa=appendix_a4_pressure_enclosure(p,x0,Tf,co,max_step),pb=appendix_a4_pressure_enclosure(p,x1,Tf,co,max_step);
  const auto zc=appendix_b_axis_data(p,m,j0,Tf,co,max_step);
  const double P0=std::max({ivabs(pa.pi0_lo,pa.pi0_hi),ivabs(pc.pi0_lo,pc.pi0_hi),ivabs(pb.pi0_lo,pb.pi0_hi)}),P1=std::max({ivabs(pa.deta_lo,pa.deta_hi),ivabs(pc.deta_lo,pc.deta_hi),ivabs(pb.deta_lo,pb.deta_hi)}),P2=std::max({ivabs(pa.detaeta_lo,pa.detaeta_hi),ivabs(pc.detaeta_lo,pc.detaeta_hi),ivabs(pb.detaeta_lo,pb.detaeta_hi)});
  const double emax=std::max(std::abs(x0),std::abs(x1)),Umax=4*emax+j0,Hmax=D*emax+(1+emax*emax)*Umax,Hpmax=D+2*emax*Umax+4*(1+emax*emax);
  const double M=appendix_b_up(invP2*A*(4*Umax*Umax+8*emax*Umax+4)+invP2*(4*Hmax+Hpmax*Umax)+2*emax*P1+P2+4*A*(P0+emax*P1));
  const double zerr=appendix_b_up(std::abs(appendix_b_d(m))*pc.deta_error+4*A*std::abs(m)*pc.pi0_error);
  glo=std::min(glo,appendix_b_down(zc.normalized_Zstar-zerr-rad*M)); ghi=std::max(ghi,appendix_b_up(zc.normalized_Zstar+zerr+rad*M)); gM=std::max(gM,M);gerr=std::max(gerr,zerr);
  finite=finite&&std::isfinite(P0)&&std::isfinite(P1)&&std::isfinite(P2);penc=penc&&pa.certified_reproduction_enclosure&&pc.certified_reproduction_enclosure&&pb.certified_reproduction_enclosure;
 }
 const double mid=.5*(a+b);const auto zmid=appendix_b_axis_data(p,mid,j0,Tf,co,max_step);const double minabs=(glo>0)?glo:((ghi<0)?-ghi:0.0);
 return {a,b,mid,zmid.normalized_Zstar,gM,glo,ghi,minabs,delta_star,gerr,finite,penc,glo>0,minabs>delta_star};
}
} // namespace nsblowup
