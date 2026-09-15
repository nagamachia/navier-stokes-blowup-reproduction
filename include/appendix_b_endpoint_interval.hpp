#pragma once

#include "appendix_b_axis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBInterval { double lo{}, hi{}; };
inline double appendix_b_down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
inline double appendix_b_up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
inline AppendixBInterval appendix_b_iv(double a,double b){return {appendix_b_down(std::min(a,b)),appendix_b_up(std::max(a,b))};}
inline double appendix_b_H_value(double eta,double h,double j0){return (0.5-h)*eta+(1.0-eta*eta)*(4.0*eta+j0);}
inline double appendix_b_Hprime_value(double eta,double h,double j0){return (4.5-h)-2.0*j0*eta-12.0*eta*eta;}

struct AppendixBEndpointIntervalAudit {
 double chi_threshold{},H_threshold{},axial_eta_lo{},axial_eta_hi{},axial_width{};
 double min_abs_H_on_azimuthal_boxes{};
 std::size_t azimuthal_boxes{},axial_boxes{};
 bool covers_eta_interval{},azimuthal_chi_certified{},axial_localized{};
};

inline AppendixBEndpointIntervalAudit appendix_b_endpoint_interval_audit(
 double h,double j0,double sigma,double eta0,double chi_threshold=0.99,
 int max_depth=80,double min_width=1e-12){
 if(!(h>0&&j0>0&&sigma>0&&chi_threshold>0&&chi_threshold<1)) throw std::invalid_argument("invalid Appendix B endpoint interval parameters");
 const double Hthr=sigma*std::sqrt(chi_threshold/(1.0-chi_threshold));
 const double disc=4*j0*j0+48*(4.5-h);
 const double rminus=(-2*j0-std::sqrt(disc))/24.0;
 const double rplus=(-2*j0+std::sqrt(disc))/24.0;
 if(!(rminus<eta0&&eta0<rplus&&appendix_b_Hprime_value(eta0,h,j0)>0)) throw std::runtime_error("central H branch not monotone");
 auto root=[&](double target,double lo,double hi){
   if(!(appendix_b_H_value(lo,h,j0)<=target&&appendix_b_H_value(hi,h,j0)>=target)) throw std::runtime_error("H threshold not bracketed");
   for(int k=0;k<max_depth && hi-lo>min_width;++k){const double m=0.5*(lo+hi); if(appendix_b_H_value(m,h,j0)<target) lo=m; else hi=m;}
   return AppendixBInterval{lo,hi};
 };
 const auto L=root(-Hthr,rminus,eta0), R=root(Hthr,eta0,rplus);
 // Include both root brackets in the axial alternative.  The complementary
 // azimuthal pieces end/start one representable number farther outward, so the
 // strict chi>threshold implication is not lost to root-rounding equality.
 const double alo=appendix_b_down(L.lo), ahi=appendix_b_up(R.hi);
 const double left_edge=appendix_b_down(alo), right_edge=appendix_b_up(ahi);
 const double leftmin=std::min({std::abs(appendix_b_H_value(-1,h,j0)),std::abs(appendix_b_H_value(rminus,h,j0)),std::abs(appendix_b_H_value(left_edge,h,j0))});
 const double rightmin=std::min({std::abs(appendix_b_H_value(rplus,h,j0)),std::abs(appendix_b_H_value(1,h,j0)),std::abs(appendix_b_H_value(right_edge,h,j0))});
 const double minabs=std::min(leftmin,rightmin);
 const bool localized=alo<=eta0&&eta0<=ahi&&(ahi-alo)<1e-4;
 const bool strict=minabs>Hthr;
 return {chi_threshold,Hthr,alo,ahi,ahi-alo,minabs,4,1,true,strict,localized};
}

} // namespace nsblowup
