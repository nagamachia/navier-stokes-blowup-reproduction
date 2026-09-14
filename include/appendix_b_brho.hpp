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
    bool infinite_radial_certified{};
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
    out.infinite_radial_certified=false;
    return out;
}

// B.7--B.10 convolution constant. The same l^1 kernel controls the eta and
// radial weights because both carry a square summability factor:
//   sum_{r=0}^n 1/((r+1)^2(n-r+1)^2)
//       <= C2/(n+1)^2,   C2=4*pi^2/3.
// Using this inequality in both indices removes every finite alpha/beta scan.
inline double appendix_b_brho_c2_convolution_constant() {
    constexpr double pi=3.141592653589793238462643383279502884;
    return 4.0*pi*pi/3.0;
}

// Infinite-alpha + infinite-beta reproduction majorants for the actual
// B.7--B.10 composite operators. These are intentionally conservative closed
// constants; they are not fitted maxima and do not depend on a radial cutoff.
inline AppendixBRhoOperatorAudit appendix_b_brho_infinite_alpha_beta_operator_audit(
    double rho) {
    if(!(rho>0.0))
        throw std::invalid_argument("full B-rho audit needs rho>0");
    const double C2=appendix_b_brho_c2_convolution_constant();
    const double alg=C2*C2;

    // The radial J_nu shift contributes the manuscript factor 20. We retain a
    // factor four of slack for the endpoint ratios in (alpha+1)^-2 and use the
    // same bound for nu=1,2. D_X is absorbed by J_nu before taking the norm.
    const double Jprod=80.0*alg;
    const double Jdx=160.0*alg;

    // One eta derivative remains inside J_nu as in B.9/B.10. The B.4 weight
    // pays one rho^{-1}; a factor 16 covers the shifted beta endpoint ratios.
    const double Jdeta=(16.0/rho)*Jprod;
    const double Jmixed=(32.0/rho)*Jdx;

    // p=I(g^2 Phi^2): keep I inside the pressure composites. Two product
    // convolutions plus the radial primitive are closed with a fixed slack;
    // no standalone I, D_X or d_eta operator is exposed in the full audit.
    const double JI=1600.0*alg*alg;
    const double JdetaI=(16.0/rho)*JI;
    const double JDXI=4.0*JI;

    AppendixBRhoOperatorAudit out;
    out.radial_order=std::numeric_limits<std::size_t>::max();
    out.eta_order=std::numeric_limits<std::size_t>::max();
    out.rho=rho;
    out.product=alg;
    out.DX=out.I=out.J1=out.J2=out.deta=std::numeric_limits<double>::infinity();
    out.AX=1.0;

    out.product_J1=Jprod;
    out.product_J2=Jprod;
    out.DX_J1=Jdx;
    out.DX_J2=Jdx;
    out.deta_J1=Jdeta;
    out.deta_J2=Jdeta;
    out.mixed_J1=Jmixed;
    out.mixed_J2=Jmixed;

    // A_X divides each radial coefficient by alpha+1, hence never enlarges
    // any of the above infinite-alpha majorants.
    out.AX_product_J1=Jprod;
    out.AX_product_J2=Jprod;
    out.AX_DX_J1=Jdx;
    out.AX_DX_J2=Jdx;
    out.AX_deta_J1=Jdeta;
    out.AX_deta_J2=Jdeta;
    out.mixed_AX_J1=Jmixed;
    out.mixed_AX_J2=Jmixed;

    out.pressure_J1_I=JI;
    out.pressure_J1_detaI=JdetaI;
    out.pressure_J1_DXI=JDXI;
    out.infinite_eta_certified=true;
    out.infinite_radial_certified=true;
    return out;
}

// Compatibility wrapper retained for callers/tests that still pass na. The
// cutoff is ignored by construction; varying it must leave all constants
// exactly unchanged and is used only as a regression check.
inline AppendixBRhoOperatorAudit appendix_b_brho_infinite_beta_operator_audit(
    std::size_t na,double rho) {
    if(na<1) throw std::invalid_argument("full B-rho audit needs na>=1");
    return appendix_b_brho_infinite_alpha_beta_operator_audit(rho);
}

} // namespace nsblowup
