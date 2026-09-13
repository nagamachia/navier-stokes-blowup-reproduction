#pragma once

#include "appendix_a7_heat_compensation.hpp"
#include "appendix_a7_log_discrepancy.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct AppendixA7LogCompensationAudit {
    double largest_target_log{};
    double smallest_target_log{};
    double target_log_spread{};
    Vector dominant_scaled_target{};
    Vector dominant_scaled_coefficients{};
    Vector dominant_scaled_linear_residual{};
    double quadratic_to_dominant_log_ratio{};
    double quadratic_log_magnitude{};
    bool all_components_resolvable_in_double{};
};

// Audit the ordered-regime correction without pretending that ordinary double
// precision can resolve every component simultaneously. The largest target is
// scaled to O(1), so the returned coefficients describe the dominant linear
// direction. Tiny target components may underflow in this representation and
// are therefore reported through target_log_spread rather than declared solved.
inline AppendixA7LogCompensationAudit appendix_a7_audit_log_compensation(
    double eta, double lambda,
    const AppendixA7SignedLogMoment& target_Cp,
    const AppendixA7SignedLogMoment& target_S,
    const AppendixA7SignedLogMoment& target_I,
    std::size_t panels = 4096) {
    const std::array<AppendixA7SignedLogMoment,3> target{
        target_Cp,target_S,target_I};
    for (const auto& q : target) {
        if (!(q.sign == -1 || q.sign == 1) || !std::isfinite(q.log_abs))
            throw std::invalid_argument("A.7 signed-log targets must be finite and nonzero");
    }

    const double Lmax = std::max({target[0].log_abs,target[1].log_abs,target[2].log_abs});
    const double Lmin = std::min({target[0].log_abs,target[1].log_abs,target[2].log_abs});
    Vector rhs(3,0.0);
    for (std::size_t i=0;i<3;++i)
        rhs[i]=static_cast<double>(target[i].sign)*std::exp(target[i].log_abs-Lmax);

    AppendixA7HeatCompensation solver(eta,lambda,panels);
    Matrix B=solver.jacobian_at_zero();

    Vector row_scale(3,0.0);
    Matrix Bs=B;
    Vector rs=rhs;
    for(std::size_t i=0;i<3;++i){
        for(double v:B[i]) row_scale[i]=std::max(row_scale[i],std::abs(v));
        if(!(row_scale[i]>0.0)) throw std::runtime_error("zero row in A.7 log Jacobian");
        for(double& v:Bs[i]) v/=row_scale[i];
        rs[i]/=row_scale[i];
    }
    const Vector c_scaled=solve_linear(Bs,rs);

    Vector residual(3,0.0);
    for(std::size_t i=0;i<3;++i){
        double v=-rhs[i];
        for(std::size_t j=0;j<3;++j) v+=B[i][j]*c_scaled[j];
        residual[i]=v;
    }

    const Vector q=solver.quadratic_remainder(c_scaled,c_scaled);
    const double qnorm=norm_inf(q);
    const double rhsnorm=std::max(1e-300,norm_inf(rhs));
    const double log_q_ratio = qnorm>0.0
        ? Lmax + std::log(qnorm/rhsnorm)
        : -INFINITY;

    AppendixA7LogCompensationAudit out;
    out.largest_target_log=Lmax;
    out.smallest_target_log=Lmin;
    out.target_log_spread=Lmax-Lmin;
    out.dominant_scaled_target=rhs;
    out.dominant_scaled_coefficients=c_scaled;
    out.dominant_scaled_linear_residual=residual;
    out.quadratic_to_dominant_log_ratio=log_q_ratio;
    out.quadratic_log_magnitude = qnorm>0.0
        ? 2.0*Lmax+std::log(qnorm)
        : -INFINITY;

    // exp(-~36) is already close to the practical component-resolution limit
    // once matrix conditioning and quadrature errors are included. This flag
    // prevents an ordered-regime audit from being mislabeled as an exact solve.
    out.all_components_resolvable_in_double = out.target_log_spread < 30.0;
    return out;
}

inline AppendixA7LogCompensationAudit appendix_a7_audit_ordered_heat_compensation(
    const OuterProfileParameters& p, double eta,
    double Tf=20.0,double co=0.05,std::size_t panels=4096) {
    const auto q=appendix_a7_leading_heat_log_discrepancy(p,eta,Tf,co,1024);
    return appendix_a7_audit_log_compensation(
        eta,p.lambda,q.patch_target_Cp,q.patch_target_S,q.patch_target_I,panels);
}

} // namespace nsblowup
