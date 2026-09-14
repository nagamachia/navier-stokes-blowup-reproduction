#pragma once

#include "appendix_b_axis.hpp"
#include "appendix_b_remainder_budget.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nsblowup {

struct AppendixBAnalyticMultiplierAudit {
    AppendixBMultiplierNorms norms;
    double rho{};
    double outer_radius{};
    double nearest_chi_zeta_pole_distance{};
    double min_H_minus_isigma{};
    double min_H_plus_isigma{};
    double qmin_H2_plus_sigma2{};
    double scalar_cauchy_factor{};
    double log_C_for_g{};
    double sup_chi{};
    double sup_zeta{};
    bool pole_separation_certified{};
    bool cauchy_certified{};
};

inline double appendix_b_distance_to_real_segment(std::complex<double> z) {
    const double x=z.real(), y=std::abs(z.imag());
    if(x<-1.0) return std::hypot(x+1.0,y);
    if(x>1.0) return std::hypot(x-1.0,y);
    return y;
}

inline std::vector<std::complex<double>> appendix_b_cubic_roots(
    std::array<std::complex<double>,4> c) {
    if(std::abs(c[3])==0.0) throw std::invalid_argument("cubic leading coefficient is zero");
    for(auto&x:c) x/=c[3];
    const double R=1.0+std::max({std::abs(c[0]),std::abs(c[1]),std::abs(c[2])});
    const std::complex<double> seed(-0.5,std::sqrt(3.0)/2.0);
    std::vector<std::complex<double>> r={R*std::complex<double>(1,0),R*seed,R*seed*seed};
    auto P=[&](std::complex<double> z){return ((z+c[2])*z+c[1])*z+c[0];};
    for(int it=0;it<250;++it){
        double md=0.0;
        const auto old=r;
        for(int i=0;i<3;++i){
            std::complex<double> den(1,0);
            for(int j=0;j<3;++j) if(i!=j) den*=old[i]-old[j];
            if(std::abs(den)<1e-30) den+=std::complex<double>(1e-30,1e-30);
            const auto d=P(old[i])/den;
            r[i]=old[i]-d;
            md=std::max(md,std::abs(d));
        }
        if(md<1e-14) break;
    }
    for(auto z:r) if(std::abs(P(z))>1e-9)
        throw std::runtime_error("Appendix B cubic root residual too large");
    return r;
}

inline double appendix_b_cauchy_brho_factor(double rho,double R) {
    if(!(rho>0.0 && R>rho)) throw std::invalid_argument("Cauchy radius must exceed rho");
    const double q=rho/R;
    double best=1.0;
    for(std::size_t b=1;b<100000;++b){
        const double term=std::pow(q,static_cast<double>(b))*std::pow(static_cast<double>(b+1),2.0);
        best=std::max(best,term);
        const double ratio=q*std::pow(static_cast<double>(b+2)/static_cast<double>(b+1),2.0);
        if(b>8 && ratio<1.0 && term<best*1e-15) break;
    }
    return best;
}

inline double appendix_b_poly_sup_from_coeffs(const std::vector<double>& c,double R) {
    const double zmax=1.0+R;
    double s=0.0,p=1.0;
    for(double a:c){s+=std::abs(a)*p;p*=zmax;}
    return s;
}

inline double appendix_b_Hprime_real_abs_bound(double h,double j0,double R) {
    // H'(x)=4.5-h-2 j0 x-12 x^2 on x in [-1-R,1+R].
    const double lo=-1.0-R, hi=1.0+R;
    auto f=[&](double x){return 4.5-h-2.0*j0*x-12.0*x*x;};
    double m=std::max(std::abs(f(lo)),std::abs(f(hi)));
    const double xv=-j0/12.0;
    if(xv>=lo && xv<=hi) m=std::max(m,std::abs(f(xv)));
    return m;
}

inline AppendixBAnalyticMultiplierAudit appendix_b_analytic_multiplier_norms(
    double h,double j0,double sigma,double rho,double radius_fraction=0.45,double Lambda=1.0) {
    if(!(h>0.0 && h<0.5) || !(j0>0.0) || !(sigma>0.0) || !(rho>0.0) ||
       !(radius_fraction>0.0 && radius_fraction<1.0) || !(Lambda>=1.0))
        throw std::invalid_argument("invalid Appendix B analytic multiplier parameters");

    const std::array<std::complex<double>,4> base={
        std::complex<double>(j0,0),std::complex<double>(4.5-h,0),
        std::complex<double>(-j0,0),std::complex<double>(-4.0,0)};
    double pole_dist=std::numeric_limits<double>::infinity();
    for(double sgn:{-1.0,1.0}){
        auto c=base;
        c[0]-=std::complex<double>(0.0,sgn*sigma);
        for(auto z:appendix_b_cubic_roots(c))
            pole_dist=std::min(pole_dist,appendix_b_distance_to_real_segment(z));
    }
    if(!(pole_dist>0.0) || !std::isfinite(pole_dist))
        throw std::runtime_error("failed to separate chi/zeta poles from [-1,1]");

    // We only need a complex neighborhood wider than rho.  Taking R=4 rho
    // keeps the all-beta Cauchy factor moderate; pole_dist remains an independent
    // numerical cross-check.  The actual no-pole proof below uses Im H directly.
    const double R=std::min(radius_fraction*pole_dist,4.0*rho);
    if(!(R>rho)) throw std::runtime_error("analytic neighborhood radius does not exceed rho");

    // For z=x+iy and H=j0+(4.5-h)z-j0 z^2-4z^3,
    //   Im H(z)=y(H'(x)+4 y^2).
    // Thus |H(z) +/- i sigma| >= sigma-|Im H(z)| uniformly on the stadium.
    const double G=appendix_b_Hprime_real_abs_bound(h,j0,R);
    const double imH=R*(G+4.0*R*R);
    const double m=sigma-imH;
    if(!(m>0.0)) throw std::runtime_error("imaginary-part separation is insufficient for chi/zeta poles");
    const double qmin=m*m;

    const double zmax=1.0+R;
    const double Hsup=appendix_b_poly_sup_from_coeffs({j0,4.5-h,-j0,-4.0},R);
    const double Lsup=1.0+2.0*h*zmax*zmax;
    const double Linf=1.0-2.0*h*zmax*zmax;
    if(!(Linf>0.0)) throw std::runtime_error("L has a zero in the chosen complex neighborhood bound");
    const double invLsup=1.0/Linf;
    const double D=0.5-h;
    const double Wsup=appendix_b_poly_sup_from_coeffs({-3.0,-2.0*D*j0,4.0-8.0*D},R);
    const double etasup=zmax;
    const double dsup=1.0+zmax*zmax;
    const double om2eUsup=appendix_b_poly_sup_from_coeffs({1.0,-2.0*j0,-8.0},R);
    const double A=0.5+h;
    const double linear_usup=appendix_b_poly_sup_from_coeffs({A+4.0,-4.0*A*j0,-16.0*A-4.0},R);
    const double Aeta4sup=4.0*A*zmax;

    // chi=1-sigma^2/Q and H/Q=(1/2)[1/(H-i sigma)+1/(H+i sigma)].
    const double chisup=1.0+sigma*sigma/qmin;
    const double zetasup=Lsup/m;
    const double F=appendix_b_cauchy_brho_factor(rho,R);

    AppendixBMultiplierNorms M{};
    M.invL=F*invLsup;
    M.Wstar=F*Wsup;
    M.Hstar=F*Hsup;
    M.eta=F*etasup;
    M.d=F*dsup;
    M.zeta=F*zetasup;
    M.chi=F*chisup;
    M.one_minus_2etaUstar=F*om2eUsup;
    M.linear_u=F*linear_usup;
    M.Aeta4=F*Aeta4sup;

    // Choose C >= sup_Omega |phi*|.  Then sup_Omega |g|<=1 and all beta
    // derivatives of g follow from the same Cauchy factor F.
    const double logC=Lambda*(1.0+R)*zetasup;

    AppendixBAnalyticMultiplierAudit out;
    out.norms=M;
    out.rho=rho;
    out.outer_radius=R;
    out.nearest_chi_zeta_pole_distance=pole_dist;
    out.min_H_minus_isigma=m;
    out.min_H_plus_isigma=m;
    out.qmin_H2_plus_sigma2=qmin;
    out.scalar_cauchy_factor=F;
    out.log_C_for_g=logC;
    out.sup_chi=chisup;
    out.sup_zeta=zetasup;
    out.pole_separation_certified=R<pole_dist && m>0.0;
    out.cauchy_certified=R>rho && std::isfinite(F) && std::isfinite(chisup) && std::isfinite(zetasup);
    return out;
}

} // namespace nsblowup
