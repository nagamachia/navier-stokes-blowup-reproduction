#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {
struct A4Interval { double lo{},hi{}; };
inline double a4_idown(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
inline double a4_iup(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
inline A4Interval a4_i(double a,double b){return {a4_idown(a),a4_iup(b)};}
inline A4Interval operator+(A4Interval x,A4Interval y){return a4_i(x.lo+y.lo,x.hi+y.hi);}
inline A4Interval operator-(A4Interval x,A4Interval y){return a4_i(x.lo-y.hi,x.hi-y.lo);}
inline A4Interval operator*(A4Interval x,A4Interval y){double q[4]={x.lo*y.lo,x.lo*y.hi,x.hi*y.lo,x.hi*y.hi};return a4_i(*std::min_element(q,q+4),*std::max_element(q,q+4));}
inline A4Interval operator*(double a,A4Interval x){return a>=0?a4_i(a*x.lo,a*x.hi):a4_i(a*x.hi,a*x.lo);}
inline A4Interval a4_inv(A4Interval x){if(x.lo<=0&&x.hi>=0)throw std::runtime_error("interval division by zero");return a4_i(1/x.hi,1/x.lo);}
inline A4Interval a4_exp(A4Interval x){return a4_i(std::exp(x.lo),std::exp(x.hi));}
inline double a4_imaxabs(A4Interval x){return std::max(std::abs(x.lo),std::abs(x.hi));}
using A4Jet=std::array<A4Interval,6>;
inline A4Jet a4_jconst(double c){A4Jet z{};for(auto&x:z)x=a4_i(0,0);z[0]=a4_i(c,c);return z;}
inline A4Jet a4_jvar(double a,double b){auto z=a4_jconst(0);z[0]=a4_i(a,b);z[1]=a4_i(1,1);return z;}
inline A4Jet a4_jscale(double c,const A4Jet&a){A4Jet z{};for(int k=0;k<=5;++k)z[k]=c*a[k];return z;}
inline A4Jet a4_jadd(const A4Jet&a,const A4Jet&b){A4Jet z{};for(int k=0;k<=5;++k)z[k]=a[k]+b[k];return z;}
inline A4Jet a4_jsub(const A4Jet&a,const A4Jet&b){A4Jet z{};for(int k=0;k<=5;++k)z[k]=a[k]-b[k];return z;}
inline A4Jet a4_jmul(const A4Jet&a,const A4Jet&b){A4Jet z{};for(int k=0;k<=5;++k){z[k]=a4_i(0,0);for(int j=0;j<=k;++j)z[k]=z[k]+a[j]*b[k-j];}return z;}
inline A4Jet a4_jinv(const A4Jet&a){A4Jet z{};z[0]=a4_inv(a[0]);for(int k=1;k<=5;++k){A4Interval s=a4_i(0,0);for(int j=1;j<=k;++j)s=s+a[j]*z[k-j];z[k]=(-1.0)*(z[0]*s);}return z;}
inline A4Jet a4_jexp(const A4Jet&a){A4Jet z{};z[0]=a4_exp(a[0]);for(int n=1;n<=5;++n){A4Interval s=a4_i(0,0);for(int k=1;k<=n;++k)s=s+(double(k)*a[k])*z[n-k];z[n]=(1.0/double(n))*s;}return z;}
inline A4Jet a4_jdiv(const A4Jet&a,const A4Jet&b){return a4_jmul(a,a4_jinv(b));}
inline A4Jet appendix_a4_smooth_step_jet(double a,double b){auto t=a4_jvar(a,b),one=a4_jconst(1.0);auto tm=a4_jsub(one,t);auto ma=a4_jscale(-1.0,a4_jinv(a4_jmul(t,t)));auto mb=a4_jscale(-1.0,a4_jinv(a4_jmul(tm,tm)));auto ea=a4_jexp(ma),eb=a4_jexp(mb);return a4_jdiv(ea,a4_jadd(ea,eb));}
struct A4SmoothStepIntervalCertificate {std::array<double,5> derivative_sup{};std::array<double,5> claimed{{9.0,100.0,4000.0,120000.0,7000000.0}};int boxes{};double endpoint_cut{};bool outward_rounded{},certified{};};
// On [eps,1/2] we use outward-rounded interval Taylor jets. The right half
// follows from s(1-t)=1-s(t). On [0,eps], exp(-1/t^2) and each of its first
// five derivatives are exponentially flat; eps=.02 makes the elementary
// exp(-1/t^2)*polynomial(1/t) bounds vastly smaller than 1e-100, while the
// denominator's other term stays >.35.  Thus 1e-100 is a conservative endpoint cap.
inline A4SmoothStepIntervalCertificate appendix_a4_smooth_step_interval_certificate(int boxes=20000,double eps=0.02){if(boxes<1||!(eps>0&&eps<0.5))throw std::invalid_argument("invalid smooth-step interval partition");std::array<double,5>M{{1e-100,1e-100,1e-100,1e-100,1e-100}};const double a=eps,b=.5,h=(b-a)/boxes,fact[6]={1,1,2,6,24,120};for(int i=0;i<boxes;++i){double x0=a+i*h,x1=(i+1==boxes)?b:a+(i+1)*h;auto j=appendix_a4_smooth_step_jet(x0,x1);for(int k=1;k<=5;++k)M[k-1]=std::max(M[k-1],fact[k]*a4_imaxabs(j[k]));}A4SmoothStepIntervalCertificate c;c.derivative_sup=M;c.boxes=boxes;c.endpoint_cut=eps;c.outward_rounded=true;c.certified=true;for(int k=0;k<5;++k)c.certified=c.certified&&std::isfinite(M[k])&&M[k]<=c.claimed[k];return c;}
} // namespace nsblowup
