#pragma once

#include "appendix_b_axis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nsblowup {

// Outward-rounded scalar interval used only for the low-degree endpoint
// geometry.  This deliberately avoids replacing the B_rho certificate by a
// point grid: every eta in [-1,1] belongs to one of the certified boxes below.
struct AppendixBInterval { double lo{}, hi{}; };

inline double appendix_b_down(double x){ return std::nextafter(x,-std::numeric_limits<double>::infinity()); }
inline double appendix_b_up(double x){ return std::nextafter(x,std::numeric_limits<double>::infinity()); }
inline AppendixBInterval appendix_b_iv(double a,double b){ return {appendix_b_down(std::min(a,b)),appendix_b_up(std::max(a,b))}; }
inline AppendixBInterval appendix_b_iadd(AppendixBInterval a,AppendixBInterval b){ return {appendix_b_down(a.lo+b.lo),appendix_b_up(a.hi+b.hi)}; }
inline AppendixBInterval appendix_b_imul(AppendixBInterval a,AppendixBInterval b){
    const double q1=a.lo*b.lo,q2=a.lo*b.hi,q3=a.hi*b.lo,q4=a.hi*b.hi;
    return {appendix_b_down(std::min(std::min(q1,q2),std::min(q3,q4))),appendix_b_up(std::max(std::max(q1,q2),std::max(q3,q4)))};
}

inline AppendixBInterval appendix_b_H_interval(double a,double b,double h,double j0){
    const auto e=appendix_b_iv(a,b);
    const auto one=appendix_b_iv(1.0,1.0);
    const auto eta2=appendix_b_imul(e,e);
    const AppendixBInterval d{appendix_b_down(1.0-eta2.hi),appendix_b_up(1.0-eta2.lo)};
    const auto U=appendix_b_iadd(appendix_b_imul(appendix_b_iv(4.0,4.0),e),appendix_b_iv(j0,j0));
    return appendix_b_iadd(appendix_b_imul(appendix_b_iv(0.5-h,0.5-h),e),appendix_b_imul(d,U));
}

struct AppendixBEndpointIntervalAudit {
    double chi_threshold{};
    double H_threshold{};
    double axial_eta_lo{};
    double axial_eta_hi{};
    double axial_width{};
    double min_abs_H_on_azimuthal_boxes{};
    std::size_t azimuthal_boxes{};
    std::size_t axial_boxes{};
    bool covers_eta_interval{};
    bool azimuthal_chi_certified{};
    bool axial_localized{};
};

// Certify the logical split chi>chi_threshold versus a small interval around
// H*=0 using outward-rounded interval evaluation of the explicit polynomial H*.
// Boxes intersecting |H|<=H_threshold are recursively bisected and retained as
// the axial alternative.  No eta sampling is used to certify the azimuthal side.
inline AppendixBEndpointIntervalAudit appendix_b_endpoint_interval_audit(
    double h,double j0,double sigma,double eta0,double chi_threshold=0.99,
    int max_depth=48,double min_width=1e-12){
    if(!(h>0.0&&j0>0.0&&sigma>0.0&&chi_threshold>0.0&&chi_threshold<1.0))
        throw std::invalid_argument("invalid Appendix B endpoint interval parameters");
    const double Hthr=sigma*std::sqrt(chi_threshold/(1.0-chi_threshold));
    struct Box{double a,b;int depth;};
    std::vector<Box> stack{{-1.0,1.0,0}};
    double alo=std::numeric_limits<double>::infinity(),ahi=-std::numeric_limits<double>::infinity();
    double minabs=std::numeric_limits<double>::infinity();
    std::size_t naz=0,nax=0;
    while(!stack.empty()){
        const auto q=stack.back(); stack.pop_back();
        const auto H=appendix_b_H_interval(q.a,q.b,h,j0);
        if(H.lo>Hthr || H.hi<-Hthr){
            ++naz; minabs=std::min(minabs,std::min(std::abs(H.lo),std::abs(H.hi))); continue;
        }
        if(q.depth>=max_depth || q.b-q.a<=min_width){
            ++nax; alo=std::min(alo,q.a); ahi=std::max(ahi,q.b); continue;
        }
        const double m=0.5*(q.a+q.b);
        stack.push_back({q.a,m,q.depth+1}); stack.push_back({m,q.b,q.depth+1});
    }
    const bool localized=nax>0 && alo<=eta0 && eta0<=ahi && (ahi-alo)<1e-4;
    return {chi_threshold,Hthr,alo,ahi,ahi-alo,minabs,naz,nax,true,naz>0&&minabs>Hthr,localized};
}

} // namespace nsblowup
