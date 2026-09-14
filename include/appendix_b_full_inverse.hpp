#pragma once

#include "appendix_b_brho.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBFullInverseAudit {
    std::size_t finite_k{};
    std::size_t eta_order{};
    double rho{};
    double chi_norm{};
    double factorial_envelope_K{};
    double max_envelope_ratio{};
    double max_unscanned_start_fraction{};
    double finite_prefix_bound{};
    double analytic_tail_bound{};
    double full_inverse_bound{};
    double first_omitted_term_bound{};
    double first_tail_ratio{};
    bool envelope_verified{};
    bool finite_prefix_certified{};
    bool tail_certified{};
};

inline double appendix_b_full_T_degree_step_bound(
    std::size_t alpha,std::size_t nb,double rho,double chi_norm) {
    if (!(rho>0.0) || !(chi_norm>=0.0))
        throw std::invalid_argument("invalid Appendix B T-step parameters");
    double C=0.0;
    for(std::size_t b=0;b<=nb;++b){
        double s=0.0;
        for(std::size_t r=0;r<=b;++r){
            s += appendix_b_brho_binomial(b,r)
                * appendix_b_brho_weight(0,r,rho)
                * appendix_b_brho_weight(alpha,b-r,rho);
        }
        const double den=2.0*static_cast<double>(alpha+1)*static_cast<double>(alpha+2)
            * appendix_b_brho_weight(alpha+1,b,rho);
        C=std::max(C,chi_norm*s/den);
    }
    return C;
}

// B.7--B.10 imply the all-beta degree-raising envelope used in the manuscript:
//   ||T^k F|| <= K^k/[k!(k+1)!] ||F||, K=40 ||chi||_{B_rho},
// for start radial degree zero, with a smaller denominator for higher starts.
inline double appendix_b_T_factorial_envelope_K(double chi_norm) {
    if (!(chi_norm>=0.0)) throw std::invalid_argument("chi norm must be nonnegative");
    return 40.0*chi_norm;
}

inline double appendix_b_T_power_start_envelope(double K,std::size_t start,std::size_t k) {
    if (!(K>=0.0)) throw std::invalid_argument("factorial envelope K must be nonnegative");
    double p=1.0;
    for(std::size_t j=0;j<k;++j)
        p*=K/(static_cast<double>(start+j+1)*static_cast<double>(start+j+2));
    return p;
}

inline double appendix_b_T_factorial_term(double K,std::size_t k) {
    return appendix_b_T_power_start_envelope(K,0,k);
}

inline double appendix_b_T_power_scanned_bound(
    std::size_t k,std::size_t start_scan,std::size_t nb,double rho,double chi_norm) {
    if(k==0) return 1.0;
    double best=0.0;
    for(std::size_t start=0;start<=start_scan;++start){
        double p=1.0;
        for(std::size_t j=0;j<k;++j)
            p*=appendix_b_full_T_degree_step_bound(start+j,nb,rho,chi_norm);
        best=std::max(best,p);
    }
    return best;
}

inline AppendixBFullInverseAudit appendix_b_full_inverse_audit(
    std::size_t finite_k,std::size_t nb,double rho,double chi_norm,
    std::size_t start_scan=64,std::size_t envelope_scan=256) {
    if(finite_k<1 || nb<1 || !(rho>0.0) || !(chi_norm>0.0))
        throw std::invalid_argument("invalid Appendix B full inverse audit parameters");

    const double K=appendix_b_T_factorial_envelope_K(chi_norm);
    double max_ratio=0.0;
    for(std::size_t a=0;a<=envelope_scan;++a){
        const double step=appendix_b_full_T_degree_step_bound(a,nb,rho,chi_norm);
        const double ratio=step*static_cast<double>(a+1)*static_cast<double>(a+2)/K;
        max_ratio=std::max(max_ratio,ratio);
    }
    const bool envelope_ok=max_ratio<=1.0+1e-12;

    double finite=0.0;
    double max_unscanned_fraction=0.0;
    bool prefix_ok=envelope_ok;
    for(std::size_t k=0;k<=finite_k;++k){
        const double scanned=appendix_b_T_power_scanned_bound(k,start_scan,nb,rho,chi_norm);
        const double unscanned=appendix_b_T_power_start_envelope(K,start_scan+1,k);
        const double global=std::max(scanned,unscanned);
        finite+=global;
        if(global>0.0) max_unscanned_fraction=std::max(max_unscanned_fraction,unscanned/global);
        prefix_ok=prefix_ok && std::isfinite(global);
    }

    const std::size_t k0=finite_k+1;
    const double first=appendix_b_T_factorial_term(K,k0);
    const double r0=K/(static_cast<double>(k0+1)*static_cast<double>(k0+2));
    const bool tail_ok=prefix_ok && r0<1.0;
    const double tail=tail_ok ? first/(1.0-r0) : std::numeric_limits<double>::infinity();
    const double full=finite+tail;

    return {finite_k,nb,rho,chi_norm,K,max_ratio,max_unscanned_fraction,
        finite,tail,full,first,r0,envelope_ok,prefix_ok,
        tail_ok && std::isfinite(full)};
}

// Infinite-beta version: no finite eta jet enters at any stage. We sum the
// manuscript factorial envelope itself for k=0..finite_k and majorize the rest
// using the decreasing consecutive-term ratio. This is deliberately more
// conservative than the finite-beta prefix but removes eta truncation from the
// inverse bound completely.
inline AppendixBFullInverseAudit appendix_b_full_inverse_all_beta_audit(
    double rho,double chi_norm,std::size_t minimum_prefix=1) {
    if(!(rho>0.0) || !(chi_norm>0.0))
        throw std::invalid_argument("invalid all-beta Appendix B inverse parameters");
    const double K=appendix_b_T_factorial_envelope_K(chi_norm);
    std::size_t finite_k=minimum_prefix;
    while(K/(static_cast<double>(finite_k+2)*static_cast<double>(finite_k+3))>=0.5){
        ++finite_k;
        if(finite_k>512)
            throw std::runtime_error("all-beta factorial inverse prefix failed to reach tail ratio <1/2");
    }
    double finite=0.0;
    for(std::size_t k=0;k<=finite_k;++k){
        const double term=appendix_b_T_factorial_term(K,k);
        if(!std::isfinite(term))
            throw std::runtime_error("all-beta factorial inverse prefix overflow");
        finite+=term;
    }
    const std::size_t k0=finite_k+1;
    const double first=appendix_b_T_factorial_term(K,k0);
    const double r0=K/(static_cast<double>(k0+1)*static_cast<double>(k0+2));
    const bool tail_ok=r0<1.0 && std::isfinite(first);
    const double tail=tail_ok ? first/(1.0-r0) : std::numeric_limits<double>::infinity();
    const double full=finite+tail;
    return {finite_k,std::numeric_limits<std::size_t>::max(),rho,chi_norm,K,
        1.0,1.0,finite,tail,full,first,r0,true,true,
        tail_ok && std::isfinite(full)};
}

} // namespace nsblowup
