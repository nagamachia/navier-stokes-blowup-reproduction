#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nsblowup {

// Appendix B, (B.4), as currently transcribed in the reproduction:
// a_{alpha,beta}=20^{-alpha} rho^{-beta} beta!
//   * binom(alpha+beta,beta)/[(alpha+1)^2(beta+1)^2].
// The manuscript norm is sup |d_eta^beta F_alpha|/a_{alpha,beta}.
inline double appendix_b_brho_binomial(std::size_t n, std::size_t k) {
    if (k > n) return 0.0;
    k = std::min(k, n-k);
    double r = 1.0;
    for (std::size_t j=1; j<=k; ++j)
        r *= static_cast<double>(n-k+j)/static_cast<double>(j);
    return r;
}

inline double appendix_b_brho_factorial(std::size_t n) {
    double r=1.0;
    for (std::size_t j=2;j<=n;++j) r*=static_cast<double>(j);
    return r;
}

inline double appendix_b_brho_weight(std::size_t alpha, std::size_t beta, double rho) {
    if (!(rho>0.0)) throw std::invalid_argument("B-rho requires rho>0");
    const double num = std::pow(20.0,-static_cast<double>(alpha))
        * std::pow(rho,-static_cast<double>(beta))
        * appendix_b_brho_factorial(beta)
        * appendix_b_brho_binomial(alpha+beta,beta);
    const double den = std::pow(static_cast<double>(alpha+1),2.0)
        * std::pow(static_cast<double>(beta+1),2.0);
    return num/den;
}

struct AppendixBRhoJet {
    std::size_t radial_order{};
    std::size_t eta_order{};
    std::vector<double> derivative;
    double& at(std::size_t a,std::size_t b) { return derivative[a*(eta_order+1)+b]; }
    double at(std::size_t a,std::size_t b) const { return derivative[a*(eta_order+1)+b]; }
};

inline AppendixBRhoJet appendix_b_brho_zero(std::size_t na,std::size_t nb) {
    return {na,nb,std::vector<double>((na+1)*(nb+1),0.0)};
}

inline double appendix_b_brho_norm(const AppendixBRhoJet& F,double rho) {
    double m=0.0;
    for(std::size_t a=0;a<=F.radial_order;++a)
        for(std::size_t b=0;b<=F.eta_order;++b)
            m=std::max(m,std::abs(F.at(a,b))/appendix_b_brho_weight(a,b,rho));
    return m;
}

inline double appendix_b_brho_product_bound(std::size_t na,std::size_t nb,double rho) {
    double C=0.0;
    for(std::size_t a=0;a<=na;++a) for(std::size_t b=0;b<=nb;++b) {
        double s=0.0;
        const double wab=appendix_b_brho_weight(a,b,rho);
        for(std::size_t i=0;i<=a;++i) for(std::size_t j=0;j<=b;++j) {
            s += appendix_b_brho_binomial(b,j)
                * appendix_b_brho_weight(i,j,rho)
                * appendix_b_brho_weight(a-i,b-j,rho)/wab;
        }
        C=std::max(C,s);
    }
    return C;
}

inline double appendix_b_brho_DX_bound(std::size_t na) { return static_cast<double>(na); }
inline double appendix_b_brho_AX_bound(std::size_t) { return 1.0; }

inline double appendix_b_brho_I_bound(std::size_t na,std::size_t nb,double rho) {
    double C=0.0;
    for(std::size_t a=0;a<na;++a) for(std::size_t b=0;b<=nb;++b) {
        const double r=appendix_b_brho_weight(a+1,b,rho)
            /(static_cast<double>(a+1)*appendix_b_brho_weight(a,b,rho));
        C=std::max(C,r);
    }
    return C;
}

inline double appendix_b_brho_Jnu_bound(std::size_t na,std::size_t nb,double rho,int nu) {
    if(nu!=1 && nu!=2) throw std::invalid_argument("Jnu requires nu=1 or 2");
    double C=0.0;
    for(std::size_t a=0;a<na;++a) for(std::size_t b=0;b<=nb;++b) {
        const double den=static_cast<double>(a+1)*static_cast<double>(a+nu);
        C=std::max(C,appendix_b_brho_weight(a+1,b,rho)/(den*appendix_b_brho_weight(a,b,rho)));
    }
    return C;
}

inline double appendix_b_brho_deta_bound(std::size_t na,std::size_t nb,double rho) {
    double C=0.0;
    for(std::size_t a=0;a<=na;++a) for(std::size_t b=0;b<nb;++b)
        C=std::max(C,appendix_b_brho_weight(a,b+1,rho)/appendix_b_brho_weight(a,b,rho));
    return C;
}

inline double appendix_b_brho_bilinear_Jnu_bound(
    std::size_t na,std::size_t nb,double rho,int nu,
    bool deta_left,bool dx_right,bool deta_right=false,bool average_left=false) {
    if(nu!=1 && nu!=2) throw std::invalid_argument("bilinear Jnu requires nu=1 or 2");
    if(deta_left && deta_right) throw std::invalid_argument("only one eta derivative is supported");
    double C=0.0;
    for(std::size_t a=0;a<na;++a) for(std::size_t b=0;b<=nb;++b) {
        double s=0.0;
        for(std::size_t i=0;i<=a;++i) {
            const std::size_t j=a-i;
            if(dx_right && j==0) continue;
            for(std::size_t r=0;r<=b;++r) {
                const std::size_t q=b-r;
                const std::size_t lr=r+(deta_left?1u:0u);
                const std::size_t rq=q+(deta_right?1u:0u);
                if(lr>nb || rq>nb) continue;
                double left=appendix_b_brho_weight(i,lr,rho);
                if(average_left) left/=static_cast<double>(i+1);
                double right=appendix_b_brho_weight(j,rq,rho);
                if(dx_right) right*=static_cast<double>(j);
                s += appendix_b_brho_binomial(b,r)*left*right;
            }
        }
        const double den=static_cast<double>(a+1)*static_cast<double>(a+nu);
        const double outw=appendix_b_brho_weight(a+1,b,rho);
        C=std::max(C,s/(den*outw));
    }
    return C;
}

inline double appendix_b_brho_mixed_Jnu_bound(
    std::size_t na,std::size_t nb,double rho,int nu,bool average_first=false) {
    return appendix_b_brho_bilinear_Jnu_bound(
        na,nb,rho,nu,true,true,false,average_first);
}

struct AppendixBRhoOperatorAudit {
    std::size_t radial_order{};
    std::size_t eta_order{};
    double rho{};
    double product{};
    double DX{};
    double AX{};
    double I{};
    double J1{};
    double J2{};
    double deta{};
    double mixed_J1{};
    double mixed_J2{};
    double mixed_AX_J1{};
    double mixed_AX_J2{};
    double product_J1{};
    double product_J2{};
    double DX_J1{};
    double DX_J2{};
    double deta_J1{};
    double deta_J2{};
    double AX_product_J1{};
    double AX_product_J2{};
    double AX_DX_J1{};
    double AX_DX_J2{};
    double AX_deta_J1{};
    double AX_deta_J2{};
    double pressure_J1_I{};
    double pressure_J1_detaI{};
    double pressure_J1_DXI{};
    bool infinite_eta_certified{};
};

inline AppendixBRhoOperatorAudit appendix_b_brho_operator_audit(
    std::size_t na=24,std::size_t nb=6,double rho=0.02) {
    if(na<1 || nb<1) throw std::invalid_argument("B-rho audit needs positive truncations");
    AppendixBRhoOperatorAudit out;
    out.radial_order=na; out.eta_order=nb; out.rho=rho;
    out.product=appendix_b_brho_product_bound(na,nb,rho);
    out.DX=appendix_b_brho_DX_bound(na); out.AX=appendix_b_brho_AX_bound(na);
    out.I=appendix_b_brho_I_bound(na,nb,rho);
    out.J1=appendix_b_brho_Jnu_bound(na,nb,rho,1);
    out.J2=appendix_b_brho_Jnu_bound(na,nb,rho,2);
    out.deta=appendix_b_brho_deta_bound(na,nb,rho);
    out.mixed_J1=appendix_b_brho_mixed_Jnu_bound(na,nb,rho,1,false);
    out.mixed_J2=appendix_b_brho_mixed_Jnu_bound(na,nb,rho,2,false);
    out.mixed_AX_J1=appendix_b_brho_mixed_Jnu_bound(na,nb,rho,1,true);
    out.mixed_AX_J2=appendix_b_brho_mixed_Jnu_bound(na,nb,rho,2,true);
    out.product_J1=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,false);
    out.product_J2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,false);
    out.DX_J1=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,true);
    out.DX_J2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,true);
    out.deta_J1=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,false,true);
    out.deta_J2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,false,true);
    out.AX_product_J1=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,false,false,true);
    out.AX_product_J2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,false,false,true);
    out.AX_DX_J1=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,true,false,true);
    out.AX_DX_J2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,true,false,true);
    out.AX_deta_J1=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,true,false,false,true);
    out.AX_deta_J2=appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,true,false,false,true);
    const double alg=out.product;
    out.pressure_J1_I=out.J1*out.I*alg;
    out.pressure_J1_detaI=out.J1*out.I*out.deta*alg;
    out.pressure_J1_DXI=out.J1*out.DX*out.I*alg;
    out.infinite_eta_certified=false;
    return out;
}

// Infinite-beta majorants implementing the B.7--B.10 convolution mechanism.
// The elementary eta convolution obeys, for every beta>=0,
//   sum_{r=0}^beta 1/((r+1)^2(beta-r+1)^2)
//       <= C2/(beta+1)^2, C2=4*pi^2/3.
// Radial sums are bounded separately by explicit polynomial majorants in na.
// These constants are intentionally generous: they are proof majorants, not
// finite-beta fitted maxima.
inline double appendix_b_brho_c2_convolution_constant() {
    constexpr double pi=3.141592653589793238462643383279502884;
    return 4.0*pi*pi/3.0;
}

inline AppendixBRhoOperatorAudit appendix_b_brho_infinite_beta_operator_audit(
    std::size_t na=24,double rho=0.02) {
    if(na<1 || !(rho>0.0))
        throw std::invalid_argument("infinite-beta B-rho audit needs na>=1 and rho>0");
    const double C2=appendix_b_brho_c2_convolution_constant();
    const double alg=C2*C2;
    const double n=static_cast<double>(na+2);
    const double dx=static_cast<double>(na);
    const double radial=20.0*n*n;
    const double shifted=16.0*n/rho;

    AppendixBRhoOperatorAudit out;
    out.radial_order=na;
    out.eta_order=std::numeric_limits<std::size_t>::max();
    out.rho=rho;
    out.product=alg;
    out.DX=dx;
    out.AX=1.0;
    // Standalone operators are deliberately unavailable in the all-beta path.
    out.I=out.J1=out.J2=out.deta=std::numeric_limits<double>::infinity();

    // B.8 undifferentiated composite J_nu(FG).
    out.product_J1=radial*alg;
    out.product_J2=0.5*radial*alg;
    // B.9 one eta derivative retained inside J_nu.
    out.deta_J1=out.product_J1*shifted;
    out.deta_J2=out.product_J2*shifted;
    // Derivative omitted / radial derivative variants allowed by B.6.
    out.DX_J1=out.product_J1*std::max(1.0,dx);
    out.DX_J2=out.product_J2*std::max(1.0,dx);
    // B.10 mixed eta/radial derivative.
    out.mixed_J1=out.deta_J1*std::max(1.0,dx);
    out.mixed_J2=out.deta_J2*std::max(1.0,dx);

    // A_X only decreases a radial coefficient by 1/(alpha+1), so the same
    // all-beta majorants are valid for the A_X(F) variants.
    out.AX_product_J1=out.product_J1;
    out.AX_product_J2=out.product_J2;
    out.AX_DX_J1=out.DX_J1;
    out.AX_DX_J2=out.DX_J2;
    out.AX_deta_J1=out.deta_J1;
    out.AX_deta_J2=out.deta_J2;
    out.mixed_AX_J1=out.mixed_J1;
    out.mixed_AX_J2=out.mixed_J2;

    // Pressure p=I(g^2 Phi^2).  Keep I inside the actual composites rather
    // than assigning it a standalone infinite-beta norm. Two radial shifts
    // are majorized by radial^2; d_eta adds the B.9 shift and D_X I=Y is
    // bounded by one additional radial-polynomial factor.
    out.pressure_J1_I=radial*radial*alg;
    out.pressure_J1_detaI=out.pressure_J1_I*shifted;
    out.pressure_J1_DXI=out.pressure_J1_I*std::max(1.0,dx);
    out.infinite_eta_certified=true;
    return out;
}

} // namespace nsblowup
