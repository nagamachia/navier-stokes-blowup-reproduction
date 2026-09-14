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
    double finite_prefix_bound{};
    double analytic_tail_bound{};
    double full_inverse_bound{};
    double first_omitted_term_bound{};
    double first_tail_ratio{};
    bool envelope_verified{};
    bool tail_certified{};
};

// Exact finite-beta B.4 one-step norm ratio for
// T=(1/2) J_2(chi .) acting on radial degree alpha.
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

// B.4 weight algebra gives the radial factorial envelope
//   ||T F_alpha|| <= K/((alpha+1)(alpha+2)) ||F_alpha||,
// with K=40 ||chi||_{B_rho}.  Consequently
//   ||T^k|| <= K^k/(k!(k+1)!)
// for the worst starting degree alpha=0.  This is the analytic radial tail
// mechanism used in the manuscript; unlike a geometric estimate it does not
// require ||T||<1.
inline double appendix_b_T_factorial_envelope_K(double chi_norm) {
    if (!(chi_norm>=0.0)) throw std::invalid_argument("chi norm must be nonnegative");
    return 40.0*chi_norm;
}

inline double appendix_b_T_factorial_term(double K,std::size_t k) {
    if (!(K>=0.0)) throw std::invalid_argument("factorial envelope K must be nonnegative");
    double term=1.0;
    for(std::size_t j=0;j<k;++j)
        term*=K/(static_cast<double>(j+1)*static_cast<double>(j+2));
    return term;
}

// Exact finite prefix, using the B.4 step constants.  We scan starting radial
// degrees because an operator norm may begin at any alpha.  The degree raising
// makes this a finite calculation for each fixed k.
inline double appendix_b_T_power_finite_bound(
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
    double finite=0.0;
    for(std::size_t k=0;k<=finite_k;++k)
        finite+=appendix_b_T_power_finite_bound(k,start_scan,nb,rho,chi_norm);

    // Verify the algebraic K/((alpha+1)(alpha+2)) envelope over a long finite
    // window as a regression guard.  The analytic bound itself is the formula
    // above; this scan detects implementation drift in the B.4 weights.
    double max_ratio=0.0;
    for(std::size_t a=0;a<=envelope_scan;++a){
        const double step=appendix_b_full_T_degree_step_bound(a,nb,rho,chi_norm);
        const double ratio=step*static_cast<double>(a+1)*static_cast<double>(a+2)/K;
        max_ratio=std::max(max_ratio,ratio);
    }
    const bool envelope_ok=max_ratio<=1.0+1e-12;

    const std::size_t k0=finite_k+1;
    const double first=appendix_b_T_factorial_term(K,k0);
    const double r0=K/(static_cast<double>(k0+1)*static_cast<double>(k0+2));
    // Ratios decrease monotonically with k.  Once r0<1, the omitted factorial
    // series is dominated by a geometric series with first term `first`.
    const bool tail_ok=envelope_ok && r0<1.0;
    const double tail=tail_ok ? first/(1.0-r0) : std::numeric_limits<double>::infinity();
    const double full=finite+tail;

    return {finite_k,nb,rho,chi_norm,K,max_ratio,finite,tail,full,first,r0,
        envelope_ok,tail_ok && std::isfinite(full)};
}

} // namespace nsblowup
