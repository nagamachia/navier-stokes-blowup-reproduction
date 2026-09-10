#include "outer_profile_stages.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {
bool close(double a, double b, double rel = 1e-10) {
    const double s = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a-b) <= rel*s;
}
}

int main() {
    using namespace nsblowup;
    if (!close(appendix_a_smooth_step(-1.0),0.0) || !close(appendix_a_smooth_step(2.0),1.0)) return 1;
    if (!close(appendix_a_smooth_step(0.5),0.5,1e-14)) return 1;
    if (!(appendix_a_smooth_step_derivative(0.5)>0.0)) return 1;

    constexpr double P=7.0, eta=0.4;
    if (!close(outer_reference_U(eta),1.6)) return 1;
    if (!close(outer_reference_E(1.0,eta,P),P/(1.0+eta*eta))) return 1;

    OuterProfileParameters params;
    params.validate();

    if (!close(outer_l_first_transition(0.0),0.6)) return 1;
    if (!close(outer_l_first_transition(1.0),0.0)) return 1;

    constexpr double Md=2.0;
    if (!close(outer_axial_k(0.0,Md),4.0)) return 1;
    const double Td=std::exp(Md)+10.0;
    if (!close(outer_axial_k(Td-10.0,Md),0.0)) return 1;
    if (!close(outer_axial_U(Td,eta,Md),0.0)) return 1;

    constexpr double lambda=0.02;
    if (!close(outer_l_intermediate_entry(0.0,lambda),0.0)) return 1;
    if (!close(outer_l_intermediate_entry(1.0,lambda),-lambda)) return 1;

    const auto p=reserved_patches_a9(lambda);
    if (!(p.profile_a < p.profile_b && p.profile_b <= p.heat_a && p.heat_b < p.positive_a && p.positive_b < p.mean_a && p.mean_a < p.mean_b)) return 1;

    constexpr double Tf=10.0, loge=2.0;
    if (!close(outer_theta_f(0.0,Tf),1.0)) return 1;
    if (!close(outer_theta_f(Tf,Tf),0.0)) return 1;
    const double final_logE=log_outer_interpolated_E(Tf,eta,loge,lambda,Tf);
    if (!close(final_logE,loge-(0.5+lambda)*Tf-std::log(2.0))) return 1;

    constexpr double h=1e-3, rho=1e-4;
    if (!(outer_terminal_f(0.0,rho)>0.0 && close(outer_terminal_f(3.0,rho),1.0))) return 1;
    for (double y : {0.0,0.5,1.0,1.5,2.0,2.5,3.0}) {
        if (!(outer_terminal_l(y,h,rho) >= -h-1e-14)) return 1;
    }
    const double Qp=outer_terminal_Qp(h,rho,240);
    if (!(Qp>0.0 && std::isfinite(Qp))) return 1;

    std::cout << "Appendix A.2 outer stage checks passed\n";
    return 0;
}
