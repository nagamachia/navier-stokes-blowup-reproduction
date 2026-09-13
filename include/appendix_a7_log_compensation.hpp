#pragma once

#include "appendix_a7_heat_compensation.hpp"
#include "appendix_a7_log_discrepancy.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixA7SignedLogCoefficient {
    int sign{};
    double log_abs{-std::numeric_limits<double>::infinity()};
};

struct AppendixA7LogCompensationSolution {
    std::array<AppendixA7SignedLogCoefficient, 3> coefficients{};
    double common_log_scale{};
    Vector scaled_coefficients{};
    Vector scaled_target{};
    Vector scaled_linear_residual{};
    double quadratic_to_linear_log_bound{};
};

inline AppendixA7LogCompensationSolution appendix_a7_solve_log_compensation(
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

    const double L = std::max({target[0].log_abs,target[1].log_abs,target[2].log_abs});
    Vector rhs(3,0.0);
    for (std::size_t i=0;i<3;++i)
        rhs[i]=static_cast<double>(target[i].sign)*std::exp(target[i].log_abs-L);

    AppendixA7HeatCompensation solver(eta,lambda,panels);
    Matrix B=solver.jacobian_at_zero();

    // Row scaling keeps the three moment equations comparably conditioned.
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

    AppendixA7LogCompensationSolution out;
    out.common_log_scale=L;
    out.scaled_coefficients=c_scaled;
    out.scaled_target=rhs;
    out.scaled_linear_residual=residual;
    for(std::size_t j=0;j<3;++j){
        if(c_scaled[j]==0.0){
            out.coefficients[j]={0,-std::numeric_limits<double>::infinity()};
        }else{
            out.coefficients[j]={c_scaled[j]>0.0?1:-1,
                                 L+std::log(std::abs(c_scaled[j]))};
        }
    }

    // Q(c,c) is homogeneous of degree two while Bc is degree one.  In the
    // common scaling c=exp(L)c_scaled, their ratio carries one extra exp(L).
    // This conservative log bound is enough to certify that the exact
    // quadratic branch is indistinguishable from its linearization in the
    // ordered regime without ever underflowing the physical coefficients.
    const Vector q=solver.quadratic_remainder(c_scaled,c_scaled);
    double qnorm=norm_inf(q);
    double bcnorm=std::max(1e-300,norm_inf(rhs));
    out.quadratic_to_linear_log_bound =
        L + (qnorm>0.0 ? std::log(qnorm/bcnorm)
                       : -std::numeric_limits<double>::infinity());
    return out;
}

inline AppendixA7LogCompensationSolution appendix_a7_solve_ordered_heat_compensation(
    const OuterProfileParameters& p, double eta,
    double Tf=20.0,double co=0.05,std::size_t panels=4096) {
    const auto q=appendix_a7_leading_heat_log_discrepancy(p,eta,Tf,co,1024);
    return appendix_a7_solve_log_compensation(
        eta,p.lambda,q.patch_target_Cp,q.patch_target_S,q.patch_target_I,panels);
}

} // namespace nsblowup
