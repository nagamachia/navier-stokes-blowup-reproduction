#pragma once

#include "appendix_a3_closure.hpp"
#include "outer_profile_stages.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace nsblowup {

struct A3NormalizedState {
    double m{};       // M/(X E)
    double j{};       // J/(X H E)
    double s{};       // S/(X E^2)
    double rI{};      // I/(X H)
    double log_E{};
};

struct A3PrePulseData {
    double eta{};
    double lambda{};
    double m_at_pulse{};
    double j_at_pulse{};
    double s_at_pulse{};
    double rI_at_pulse{};
    double pre_M_at_first_bump{};
    double pre_J_at_first_bump{};
    double log_E_at_pulse{};
};

struct A3GlobalClosure {
    double amplitude{};
    double s_infinity{}; // S(infinity)/(Xp eb^2 f^2)
    double pulse_residual_M{};
    double pulse_residual_J{};
    double rI_exterior_start{};
    double rI_target{};
    double angular_c1{};
    double angular_c2{};
    double exterior_hold_length{};
    double terminal_weight{}; // (X E^2)/(Xp eb^2 f^2) at end of terminal interval
};

namespace detail {

inline A3NormalizedState a3_rhs_ratio(const A3NormalizedState& q,
                                      double u_over_e,
                                      double l) {
    return {
        u_over_e - (0.5 + l) * q.m,
        u_over_e - (0.5 + 2.0 * l) * q.j,
        u_over_e * u_over_e - 0.5 - 2.0 * l * q.s,
        1.0 - (1.0 + l) * q.rI,
        l - 0.5};
}

inline A3NormalizedState a3_rhs(const A3NormalizedState& q,
                                double U,
                                double l) {
    const double u_over_e = U == 0.0 ? 0.0 : U * std::exp(-q.log_E);
    return a3_rhs_ratio(q,u_over_e,l);
}

inline A3NormalizedState add_scaled(const A3NormalizedState& a,
                                    const A3NormalizedState& b,
                                    double scale) {
    return {a.m + scale*b.m,
            a.j + scale*b.j,
            a.s + scale*b.s,
            a.rI + scale*b.rI,
            a.log_E + scale*b.log_E};
}

template<class UFunc, class LFunc>
A3NormalizedState integrate_stage(A3NormalizedState q,
                                  double length,
                                  UFunc&& U,
                                  LFunc&& l,
                                  double max_step = 0.02) {
    if (length <= 0.0) return q;
    const int steps = std::max(1, static_cast<int>(std::ceil(length/max_step)));
    const double h = length/static_cast<double>(steps);
    for (int n=0; n<steps; ++n) {
        const double y = n*h;
        const auto k1=a3_rhs(q,U(y),l(y));
        const auto q2=add_scaled(q,k1,0.5*h);
        const auto k2=a3_rhs(q2,U(y+0.5*h),l(y+0.5*h));
        const auto q3=add_scaled(q,k2,0.5*h);
        const auto k3=a3_rhs(q3,U(y+0.5*h),l(y+0.5*h));
        const auto q4=add_scaled(q,k3,h);
        const auto k4=a3_rhs(q4,U(y+h),l(y+h));
        q.m += h*(k1.m+2*k2.m+2*k3.m+k4.m)/6.0;
        q.j += h*(k1.j+2*k2.j+2*k3.j+k4.j)/6.0;
        q.s += h*(k1.s+2*k2.s+2*k3.s+k4.s)/6.0;
        q.rI += h*(k1.rI+2*k2.rI+2*k3.rI+k4.rI)/6.0;
        q.log_E += h*(k1.log_E+2*k2.log_E+2*k3.log_E+k4.log_E)/6.0;
    }
    return q;
}

inline double bump_derivative(double y,double center,double width=0.3) {
    const double eps=1e-5*width;
    return (appendix_a3_unit_bump(y+eps,center,width)-
            appendix_a3_unit_bump(y-eps,center,width))/(2.0*eps);
}

inline double profile_interpolation_l(double y,double eta,double lambda,double Tf) {
    const double theta_prime=-appendix_a_smooth_step_derivative(y/Tf)/Tf;
    return -lambda + theta_prime*(std::log(2.0)-std::log1p(eta*eta));
}

inline double advance_rI_constant_l(double r,double l,double length) {
    const double k=1.0+l;
    if (std::abs(k)<1e-14) return r+length;
    const double equilibrium=1.0/k;
    return equilibrium+(r-equilibrium)*std::exp(-k*length);
}

struct SWeightState { double sbar{}; double weight{}; };

inline SWeightState advance_S_constant_l(SWeightState q,double l,double length) {
    if (length<=0.0) return q;
    if (std::abs(l)<1e-14) {
        q.sbar-=0.5*q.weight*length;
        return q;
    }
    const double factor=std::exp(2.0*l*length);
    q.sbar-=q.weight*std::expm1(2.0*l*length)/(4.0*l);
    q.weight*=factor;
    return q;
}

template<class LFunc>
SWeightState integrate_S_weight(SWeightState q,double length,LFunc&& l,double max_step=0.02) {
    if (length<=0.0) return q;
    const int steps=std::max(1,static_cast<int>(std::ceil(length/max_step)));
    const double h=length/static_cast<double>(steps);
    for(int n=0;n<steps;++n){
        const double y=n*h;
        auto deriv=[&](const SWeightState& z,double yy){
            const double ll=l(yy);
            return SWeightState{-0.5*z.weight,2.0*ll*z.weight};
        };
        const auto k1=deriv(q,y);
        SWeightState z2{q.sbar+0.5*h*k1.sbar,q.weight+0.5*h*k1.weight};
        const auto k2=deriv(z2,y+0.5*h);
        SWeightState z3{q.sbar+0.5*h*k2.sbar,q.weight+0.5*h*k2.weight};
        const auto k3=deriv(z3,y+0.5*h);
        SWeightState z4{q.sbar+h*k3.sbar,q.weight+h*k3.weight};
        const auto k4=deriv(z4,y+h);
        q.sbar+=h*(k1.sbar+2*k2.sbar+2*k3.sbar+k4.sbar)/6.0;
        q.weight+=h*(k1.weight+2*k2.weight+2*k3.weight+k4.weight)/6.0;
    }
    return q;
}

template<class LFunc>
double integrate_Q(double Q,double length,double hpar,LFunc&& l,double max_step=0.02){
    if(length<=0.0) return Q;
    const int steps=std::max(1,static_cast<int>(std::ceil(length/max_step)));
    const double ds=length/static_cast<double>(steps);
    auto rhs=[&](double q,double y){const double ll=l(y); return -(1.0+ll)*q-ll-hpar;};
    for(int n=0;n<steps;++n){
        const double y=n*ds;
        const double k1=rhs(Q,y);
        const double k2=rhs(Q+0.5*ds*k1,y+0.5*ds);
        const double k3=rhs(Q+0.5*ds*k2,y+0.5*ds);
        const double k4=rhs(Q+ds*k3,y+ds);
        Q+=ds*(k1+2*k2+2*k3+k4)/6.0;
    }
    return Q;
}

} // namespace detail

inline A3PrePulseData appendix_a3_integrate_to_pulse(const OuterProfileParameters& p,
                                                     double eta,
                                                     double max_step = 0.02) {
    p.validate();
    if (!(eta >= -1.0 && eta <= 1.0))
        throw std::invalid_argument("eta must lie in [-1,1]");

    const double f=outer_shape_f(eta);
    const double logE0=p.log_Pstar+std::log(f);
    const double invE=std::exp(-logE0);

    A3NormalizedState q;
    q.m=4.0*eta*invE;
    q.j=2.5*eta*invE;
    q.s=16.0*eta*eta*invE*invE-5.0/12.0;
    q.rI=1.0/1.6;
    q.log_E=logE0;

    q=detail::integrate_stage(q,1.0,
        [=](double){return 4.0*eta;},
        [](double y){return outer_l_first_transition(y);},max_step);
    q=detail::integrate_stage(q,p.Td(),
        [&](double y){return outer_axial_U(y,eta,p.Md);},
        [](double){return 0.0;},max_step);
    q=detail::integrate_stage(q,1.0,
        [](double){return 0.0;},
        [&](double y){return outer_l_intermediate_entry(y,p.lambda);},max_step);
    q=detail::integrate_stage(q,p.Tw(),
        [](double){return 0.0;},
        [&](double){return -p.lambda;},max_step);

    const double y1=13.0/p.lambda-3.0;
    const double s1=0.5-p.lambda;
    const double s2=0.5-2.0*p.lambda;

    A3PrePulseData out;
    out.eta=eta;
    out.lambda=p.lambda;
    out.m_at_pulse=q.m;
    out.j_at_pulse=q.j;
    out.s_at_pulse=q.s;
    out.rI_at_pulse=q.rI;
    out.pre_M_at_first_bump=q.m*std::exp(-s1*y1);
    out.pre_J_at_first_bump=q.j*std::exp(-s2*y1);
    out.log_E_at_pulse=q.log_E;
    return out;
}

inline PulseMomentClosure appendix_a3_close_MJ_from_schedule(
    const OuterProfileParameters& p,double eta,double amplitude,double max_step=0.02) {
    const auto pre=appendix_a3_integrate_to_pulse(p,eta,max_step);
    return appendix_a3_close_pulse_MJ(p.lambda,pre.pre_M_at_first_bump,
                                      pre.pre_J_at_first_bump,amplitude);
}

// Directly evaluate S(infinity) in the fixed pulse-start normalization
// Xp*eb^2*f^2.  This avoids catastrophic scale growth across the h-dependent
// exterior release and turns the remainder in (A.19) into a computed quantity.
inline A3GlobalClosure appendix_a3_evaluate_global_closure(
    const OuterProfileParameters& p,double eta,double amplitude,
    double Tf=20.0,double co=0.05,double max_step=0.03) {
    p.validate();
    const double h=std::exp(p.log_h);
    if (!(h>0.0)) throw std::runtime_error("h underflowed; increase log_h for numerical evaluation");
    if (!(Tf>0.0 && co>0.0)) throw std::invalid_argument("Tf and co must be positive");

    const auto pre=appendix_a3_integrate_to_pulse(p,eta,max_step);
    const auto mj=appendix_a3_close_pulse_MJ(p.lambda,pre.pre_M_at_first_bump,
                                             pre.pre_J_at_first_bump,amplitude);

    // Pulse contribution in xi=lambda*y.  Integrating in y resolves the two
    // width-.3 end bumps without introducing a lambda-dependent quadrature grid.
    const double pulse_len=13.0/p.lambda;
    const double y1=pulse_len-3.0, y2=pulse_len-1.0;
    const int n=std::max(2,static_cast<int>(std::ceil(pulse_len/max_step/2.0))*2);
    const double dy=pulse_len/static_cast<double>(n);
    auto pulse_integrand=[&](double y){
        const double R=amplitude*appendix_a3_R0(p.lambda*y)+
                       mj.c1*appendix_a3_unit_bump(y,y1)+
                       mj.c2*appendix_a3_unit_bump(y,y2);
        return std::exp(-2.0*p.lambda*y)*(R*R-0.5);
    };
    double pulse_sum=pulse_integrand(0.0)+pulse_integrand(pulse_len);
    for(int i=1;i<n;++i) pulse_sum+=(i%2?4.0:2.0)*pulse_integrand(i*dy);
    double sbar=pre.s_at_pulse+pulse_sum*dy/3.0;
    double weight=std::exp(-26.0);

    // I/(XH) is independent of U, so propagate it separately through the pulse.
    double rI=detail::advance_rI_constant_l(pre.rI_at_pulse,-p.lambda,pulse_len);

    // A.10 profile interpolation.
    auto lint=[&](double y){return detail::profile_interpolation_l(y,eta,p.lambda,Tf);};
    detail::SWeightState sw{sbar,weight};
    sw=detail::integrate_S_weight(sw,Tf,lint,max_step);
    A3NormalizedState rq{0,0,0,rI,0};
    rq=detail::integrate_stage(rq,Tf,[](double){return 0.0;},lint,max_step);
    rI=rq.rI;

    // Uniform A.10 hold until the first A.11 bump center.
    const double uniform_len=30.0*std::log(1.0/p.lambda);
    const double first_center=uniform_len-3.0;
    rI=detail::advance_rI_constant_l(rI,-p.lambda,first_center);
    sw=detail::advance_S_constant_l(sw,-p.lambda,first_center);

    const double target_rI=1.0/(1.0-p.lambda);
    const double dI=target_rI-rI;
    const auto ip=appendix_a3_close_angular_I_pressure(p.lambda,dI,0.0);

    auto corrected_l=[&](double t){
        const double y=first_center+t;
        const double b1=appendix_a3_unit_bump(y,first_center);
        const double b2=appendix_a3_unit_bump(y,uniform_len-1.0);
        const double factor=1.0+ip.c1*b1+ip.c2*b2;
        if (!(factor>0.0)) throw std::runtime_error("A.11 relative E correction lost positivity");
        const double fp=ip.c1*detail::bump_derivative(y,first_center)+
                        ip.c2*detail::bump_derivative(y,uniform_len-1.0);
        return -p.lambda+fp/factor;
    };
    sw=detail::integrate_S_weight(sw,3.0,corrected_l,std::min(max_step,0.01));
    A3NormalizedState rq2{0,0,0,rI,0};
    rq2=detail::integrate_stage(rq2,3.0,[](double){return 0.0;},corrected_l,std::min(max_step,0.01));
    rI=rq2.rI;

    // Exterior release: -lambda -> -1, hold -1, then -1 -> -h.
    auto ldown=[&](double y){return -p.lambda+(-1.0+p.lambda)*appendix_a_smooth_step(y);};
    sw=detail::integrate_S_weight(sw,1.0,ldown,max_step);
    double Q=(p.lambda-h)/(1.0-p.lambda);
    Q=detail::integrate_Q(Q,1.0,h,ldown,max_step);

    const double hold1=4.0*std::log(1.0/h);
    sw=detail::advance_S_constant_l(sw,-1.0,hold1);
    Q+=(1.0-h)*hold1; // exact for l=-1

    auto lup=[&](double y){return -1.0+(1.0-h)*appendix_a_smooth_step(y);};
    sw=detail::integrate_S_weight(sw,1.0,lup,max_step);
    Q=detail::integrate_Q(Q,1.0,h,lup,max_step);

    const double rho_o=co*h;
    const double Qp=outer_terminal_Qp(h,rho_o,600);
    if (!(Q>Qp && Qp>0.0)) throw std::runtime_error("exterior Q hold is not well ordered");
    const double hold_h=std::log(Q/Qp)/(1.0-h);
    sw=detail::advance_S_constant_l(sw,-h,hold_h);

    auto lterm=[&](double y){return outer_terminal_l(y,h,rho_o);};
    sw=detail::integrate_S_weight(sw,3.0,lterm,std::min(max_step,0.01));

    // Beyond the terminal interval E=Epow ~ X^{-1/2-h}; add its exact tail.
    const double tail=-sw.weight/(4.0*h);
    sbar=sw.sbar+tail;

    A3GlobalClosure out;
    out.amplitude=amplitude;
    out.s_infinity=sbar;
    out.pulse_residual_M=mj.residual_M;
    out.pulse_residual_J=mj.residual_J;
    out.rI_exterior_start=rI;
    out.rI_target=target_rI;
    out.angular_c1=ip.c1;
    out.angular_c2=ip.c2;
    out.exterior_hold_length=hold_h;
    out.terminal_weight=sw.weight;
    return out;
}

inline A3GlobalClosure appendix_a3_solve_global_amplitude(
    const OuterProfileParameters& p,double eta,double Tf=20.0,double co=0.05,
    double max_step=0.03,int iterations=60) {
    auto value=[&](double a){return appendix_a3_evaluate_global_closure(p,eta,a,Tf,co,max_step);};
    double lo=0.9,hi=1.2;
    auto vlo=value(lo), vhi=value(hi);
    if (!(vlo.s_infinity<0.0 && vhi.s_infinity>0.0))
        throw std::runtime_error("global A.3 S(infinity) root not bracketed on [0.9,1.2]");
    for(int i=0;i<iterations;++i){
        const double mid=0.5*(lo+hi);
        const auto vm=value(mid);
        if(vm.s_infinity>0.0) hi=mid; else lo=mid;
    }
    return value(0.5*(lo+hi));
}

} // namespace nsblowup
