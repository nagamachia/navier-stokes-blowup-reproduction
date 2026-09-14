#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nsblowup {

// Appendix B, (B.4):
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

// Direct finite-truncation B.6 family.  The flags select
// J_nu[(d_eta F)(D_X G)], J_nu[F(D_X G)], J_nu[F d_eta G], or J_nu[FG].
// average_left replaces F by A_X(F).  Keeping J_nu inside this sum is the
// essential point: splitting d_eta and D_X into separate operator norms loses
// the alpha-beta compensation built into (B.4).
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
};

inline AppendixBRhoOperatorAudit appendix_b_brho_operator_audit(
    std::size_t na=24,std::size_t nb=6,double rho=0.02) {
    if(na<1 || nb<1) throw std::invalid_argument("B-rho audit needs positive truncations");
    return {na,nb,rho,
        appendix_b_brho_product_bound(na,nb,rho),
        appendix_b_brho_DX_bound(na), appendix_b_brho_AX_bound(na),
        appendix_b_brho_I_bound(na,nb,rho),
        appendix_b_brho_Jnu_bound(na,nb,rho,1),
        appendix_b_brho_Jnu_bound(na,nb,rho,2),
        appendix_b_brho_deta_bound(na,nb,rho),
        appendix_b_brho_mixed_Jnu_bound(na,nb,rho,1,false),
        appendix_b_brho_mixed_Jnu_bound(na,nb,rho,2,false),
        appendix_b_brho_mixed_Jnu_bound(na,nb,rho,1,true),
        appendix_b_brho_mixed_Jnu_bound(na,nb,rho,2,true),
        appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,false),
        appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,false),
        appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,true),
        appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,true),
        appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,1,false,false,true),
        appendix_b_brho_bilinear_Jnu_bound(na,nb,rho,2,false,false,true)};
}

} // namespace nsblowup
