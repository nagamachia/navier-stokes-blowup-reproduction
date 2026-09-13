#pragma once

#include "appendix_b_axis.hpp"
#include "appendix_b_remainder_budget.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace nsblowup {

// Truncated Taylor jet: c[k]=f^(k)(eta)/k!.  Algebra is exact up to order n.
struct AppendixBScalarJet {
    std::vector<double> c;
};
inline AppendixBScalarJet appendix_b_sj_const(double x,std::size_t n){
    AppendixBScalarJet a{std::vector<double>(n+1,0.0)}; a.c[0]=x; return a;
}
inline AppendixBScalarJet appendix_b_sj_var(double x,std::size_t n){
    auto a=appendix_b_sj_const(x,n); if(n) a.c[1]=1.0; return a;
}
inline AppendixBScalarJet appendix_b_sj_add(const AppendixBScalarJet&a,const AppendixBScalarJet&b){
    AppendixBScalarJet r{std::vector<double>(a.c.size(),0.0)};
    for(std::size_t k=0;k<r.c.size();++k) r.c[k]=a.c[k]+b.c[k]; return r;
}
inline AppendixBScalarJet appendix_b_sj_scale(const AppendixBScalarJet&a,double s){
    auto r=a; for(double&x:r.c)x*=s; return r;
}
inline AppendixBScalarJet appendix_b_sj_mul(const AppendixBScalarJet&a,const AppendixBScalarJet&b){
    AppendixBScalarJet r{std::vector<double>(a.c.size(),0.0)};
    for(std::size_t k=0;k<r.c.size();++k) for(std::size_t j=0;j<=k;++j) r.c[k]+=a.c[j]*b.c[k-j]; return r;
}
inline AppendixBScalarJet appendix_b_sj_inv(const AppendixBScalarJet&a){
    if(a.c.empty() || a.c[0]==0.0) throw std::invalid_argument("singular scalar jet inverse");
    AppendixBScalarJet r{std::vector<double>(a.c.size(),0.0)}; r.c[0]=1.0/a.c[0];
    for(std::size_t k=1;k<r.c.size();++k){
        double s=0.0; for(std::size_t j=1;j<=k;++j)s+=a.c[j]*r.c[k-j]; r.c[k]=-s/a.c[0];
    } return r;
}
inline AppendixBScalarJet appendix_b_sj_div(const AppendixBScalarJet&a,const AppendixBScalarJet&b){ return appendix_b_sj_mul(a,appendix_b_sj_inv(b)); }

inline double appendix_b_scalar_brho_norm(const AppendixBScalarJet& a,double rho){
    double m=0.0, fact=1.0, rp=1.0;
    for(std::size_t b=0;b<a.c.size();++b){
        if(b>0){fact*=static_cast<double>(b); rp*=rho;}
        const double deriv=std::abs(a.c[b])*fact;
        const double w=std::pow(rho,-static_cast<double>(b))*fact/std::pow(static_cast<double>(b+1),2.0);
        m=std::max(m,deriv/w);
    }
    return m;
}

struct AppendixBMultiplierSampleAudit {
    AppendixBMultiplierNorms norms;
    double rho{};
    std::size_t eta_order{};
    double sigma_star{};
    double eta0{};
    int samples{};
};

inline AppendixBMultiplierSampleAudit appendix_b_sample_multiplier_norms(
    double h,double j0,double sigma,double eta0,double rho,std::size_t nb=6,int samples=4001){
    if(!(h>0.0 && h<0.5) || !(sigma>0.0) || !(rho>0.0) || samples<101)
        throw std::invalid_argument("invalid multiplier audit parameters");
    AppendixBMultiplierNorms M{};
    auto audit=[&](double eta){
        const auto e=appendix_b_sj_var(eta,nb), one=appendix_b_sj_const(1.0,nb);
        const auto e2=appendix_b_sj_mul(e,e);
        const auto d=appendix_b_sj_add(one,appendix_b_sj_scale(e2,-1.0));
        const auto L=appendix_b_sj_add(one,appendix_b_sj_scale(e2,-2.0*h));
        const auto invL=appendix_b_sj_inv(L);
        const auto U=appendix_b_sj_add(appendix_b_sj_scale(e,4.0),appendix_b_sj_const(j0,nb));
        const auto H=appendix_b_sj_add(appendix_b_sj_scale(e,0.5-h),appendix_b_sj_mul(d,U));
        const auto W=appendix_b_sj_add(one,
            appendix_b_sj_add(appendix_b_sj_scale(d,-4.0),
                appendix_b_sj_scale(appendix_b_sj_mul(e,U),-2.0*(0.5-h))));
        const auto H2=appendix_b_sj_mul(H,H);
        const auto denom=appendix_b_sj_add(H2,appendix_b_sj_const(sigma*sigma,nb));
        const auto zeta=appendix_b_sj_scale(appendix_b_sj_div(appendix_b_sj_mul(L,H),denom),-1.0);
        const auto om2eU=appendix_b_sj_add(one,appendix_b_sj_scale(appendix_b_sj_mul(e,U),-2.0));
        const auto lin=appendix_b_sj_add(
            appendix_b_sj_scale(appendix_b_sj_add(one,appendix_b_sj_scale(appendix_b_sj_mul(e,U),-4.0)),0.5+h),
            appendix_b_sj_scale(d,4.0));
        const auto Aeta4=appendix_b_sj_scale(e,4.0*(0.5+h));
        M.invL=std::max(M.invL,appendix_b_scalar_brho_norm(invL,rho));
        M.Wstar=std::max(M.Wstar,appendix_b_scalar_brho_norm(W,rho));
        M.Hstar=std::max(M.Hstar,appendix_b_scalar_brho_norm(H,rho));
        M.eta=std::max(M.eta,appendix_b_scalar_brho_norm(e,rho));
        M.d=std::max(M.d,appendix_b_scalar_brho_norm(d,rho));
        M.zeta=std::max(M.zeta,appendix_b_scalar_brho_norm(zeta,rho));
        M.one_minus_2etaUstar=std::max(M.one_minus_2etaUstar,appendix_b_scalar_brho_norm(om2eU,rho));
        M.linear_u=std::max(M.linear_u,appendix_b_scalar_brho_norm(lin,rho));
        M.Aeta4=std::max(M.Aeta4,appendix_b_scalar_brho_norm(Aeta4,rho));
    };
    for(int i=0;i<samples;++i) audit(-1.0+2.0*static_cast<double>(i)/static_cast<double>(samples-1));
    // zeta derivatives peak on the sigma/H' scale near H*=0; resolve that scale explicitly.
    const double Hprime=std::abs((0.5-h)+4.0*(1.0-eta0*eta0)-2.0*eta0*(4.0*eta0+j0));
    const double scale=sigma/std::max(1e-14,Hprime);
    for(int k=-80;k<=80;++k){ const double x=eta0+0.125*static_cast<double>(k)*scale; if(x>=-1.0&&x<=1.0)audit(x); }
    return {M,rho,nb,sigma,eta0,samples};
}

} // namespace nsblowup
