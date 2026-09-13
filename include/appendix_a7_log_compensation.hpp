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
    double inverse_jacobian_inf_bound{};
    double quadratic_operator_inf_bound{};
    double log_contraction_factor_bound{};
    double log_fixed_point_radius_bound{};
    bool nonlinear_branch_contraction_certified{};
    bool all_components_resolvable_in_double{};
};

// Audit the ordered-regime correction without pretending that ordinary double
// precision can resolve every component simultaneously. The largest target is
// scaled to O(1), so the returned coefficients describe the dominant linear
// direction. Tiny target components may underflow in this representation and
// are therefore reported through target_log_spread rather than declared solved.
//
// For the exact nonlinear branch write c=eps*z, eps=exp(Lmax). After row
// scaling the equation is
//     Bs z + eps Qs(z,z) = d,  ||d||_inf = 1.
// Let beta >= ||Bs^{-1}||_inf and q >= ||Qs||_{inf,inf->inf}. On the ball
// ||z||<=2 beta, the zero-start fixed-point map is invariant and contractive
// whenever 4*eps*beta^2*q < 1. We store the logarithm of this quantity so the
// paper's exponentially ordered regime never underflows.
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

    const Vector q_at_c=solver.quadratic_remainder(c_scaled,c_scaled);
    const double qnorm=norm_inf(q_at_c);
    const double rhsnorm=std::max(1e-300,norm_inf(rhs));
    const double log_q_ratio = qnorm>0.0
        ? Lmax + std::log(qnorm/rhsnorm)
        : -INFINITY;

    // Compute ||Bs^{-1}||_inf from its three columns.
    Matrix inv(3,Vector(3,0.0));
    for(std::size_t k=0;k<3;++k){
        Vector e(3,0.0); e[k]=1.0;
        const Vector col=solve_linear(Bs,e);
        for(std::size_t i=0;i<3;++i) inv[i][k]=col[i];
    }
    double beta=0.0;
    for(std::size_t i=0;i<3;++i){
        double rowsum=0.0;
        for(double v:inv[i]) rowsum+=std::abs(v);
        beta=std::max(beta,rowsum);
    }

    // Bound the row-scaled bilinear map by summing absolute coefficients
    // Qs_i(e_j,e_k). This is a direct finite-dimensional operator-norm bound.
    double qbound=0.0;
    for(std::size_t i=0;i<3;++i){
        double rowsum=0.0;
        for(std::size_t j=0;j<3;++j){
            Vector ej(3,0.0); ej[j]=1.0;
            for(std::size_t k=0;k<3;++k){
                Vector ek(3,0.0); ek[k]=1.0;
                const Vector qjk=solver.quadratic_remainder(ej,ek);
                rowsum+=std::abs(qjk[i]/row_scale[i]);
            }
        }
        qbound=std::max(qbound,rowsum);
    }
    if (!(beta>0.0 && qbound>0.0 && std::isfinite(beta) && std::isfinite(qbound)))
        throw std::runtime_error("invalid A.7 contraction constants");

    const double log_contraction = Lmax + std::log(4.0) +
        2.0*std::log(beta) + std::log(qbound);
    const double log_radius = Lmax + std::log(2.0*beta);

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
    out.inverse_jacobian_inf_bound=beta;
    out.quadratic_operator_inf_bound=qbound;
    out.log_contraction_factor_bound=log_contraction;
    out.log_fixed_point_radius_bound=log_radius;
    out.nonlinear_branch_contraction_certified=log_contraction<0.0;

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
