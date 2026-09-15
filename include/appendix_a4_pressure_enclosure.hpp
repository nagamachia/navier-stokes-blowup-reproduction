#pragma once

#include "appendix_a4_pressure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct A4PressureEnclosure {
    A4PressureDatum center{};
    double pi0_error{};
    double deta_error{};
    double detaeta_error{};
    double pi0_lo{},pi0_hi{};
    double deta_lo{},deta_hi{};
    double detaeta_lo{},detaeta_hi{};
    double coarse_step{};
    double fine_step{};
    bool nested_convergence{};
    bool certified_reproduction_enclosure{};
};

inline double a4_enclosure_down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
inline double a4_enclosure_up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}

// A posteriori radial quadrature enclosure used by the numerical reproduction.
// The schedule consists of exact constant-l pieces plus RK4 transition pieces
// and composite Simpson on A.10.  Both variable-step contributions are fourth
// order globally.  For two nested steps H and H/2 we therefore use the standard
// Richardson defect |fine-coarse|/(2^4-1), enlarged by a safety factor and by a
// second H/4 comparison.  This is a validated numerical-reproduction enclosure,
// not yet a formal theorem about floating-point RK4: a future interval-ODE pass
// can replace this estimator without changing downstream APIs.
inline A4PressureEnclosure appendix_a4_pressure_enclosure(
    const OuterProfileParameters& p,double eta,double Tf=20.0,double co=0.05,
    double coarse_step=0.04,double safety=8.0){
    if(!(coarse_step>0.0&&safety>=1.0)) throw std::invalid_argument("invalid A.4 enclosure parameters");
    const auto q0=appendix_a4_pressure_datum(p,eta,Tf,co,coarse_step);
    const auto q1=appendix_a4_pressure_datum(p,eta,Tf,co,0.5*coarse_step);
    const auto q2=appendix_a4_pressure_datum(p,eta,Tf,co,0.25*coarse_step);
    auto err=[&](double a,double b,double c){
        const double e01=std::abs(b-a)/15.0;
        const double e12=std::abs(c-b)/15.0;
        return a4_enclosure_up(safety*std::max(e12,e01/16.0));
    };
    const double e0=err(q0.normalized_pi0,q1.normalized_pi0,q2.normalized_pi0);
    const double e1=err(q0.normalized_deta,q1.normalized_deta,q2.normalized_deta);
    const double e2=err(q0.normalized_detaeta,q1.normalized_detaeta,q2.normalized_detaeta);
    const bool nested = std::abs(q2.normalized_pi0-q1.normalized_pi0) <= 1.25*std::abs(q1.normalized_pi0-q0.normalized_pi0)+1e-14
                     && std::abs(q2.normalized_deta-q1.normalized_deta) <= 1.25*std::abs(q1.normalized_deta-q0.normalized_deta)+1e-14
                     && std::abs(q2.normalized_detaeta-q1.normalized_detaeta) <= 1.25*std::abs(q1.normalized_detaeta-q0.normalized_detaeta)+1e-14;
    return {q2,e0,e1,e2,
      a4_enclosure_down(q2.normalized_pi0-e0),a4_enclosure_up(q2.normalized_pi0+e0),
      a4_enclosure_down(q2.normalized_deta-e1),a4_enclosure_up(q2.normalized_deta+e1),
      a4_enclosure_down(q2.normalized_detaeta-e2),a4_enclosure_up(q2.normalized_detaeta+e2),
      coarse_step,0.25*coarse_step,nested,nested&&std::isfinite(e0)&&std::isfinite(e1)&&std::isfinite(e2)};
}

} // namespace nsblowup
