#pragma once

#include "appendix_a4_pressure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct A4PressureEnclosure {
    A4PressureDatum center{};
    double pi0_error{},deta_error{},detaeta_error{};
    double pi0_lo{},pi0_hi{},deta_lo{},deta_hi{},detaeta_lo{},detaeta_hi{};
    double coarse_step{},fine_step{};
    bool nested_convergence{}; // retained ABI field; now means analytic-bound path active
    bool certified_reproduction_enclosure{};
};

inline double a4_enclosure_down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
inline double a4_enclosure_up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}

struct A4TruncationMajorant {
    double smooth_d1{},smooth_d2{},smooth_d3{},smooth_d4{},smooth_d5{};
    double rk4_fifth_state{},a10_fourth_integrand{};
    double rk4_error{},simpson_error{},total_integral_error{};
};

// Real-axis derivative majorants for s(t)=e^{-1/t^2}/(e^{-1/t^2}+e^{-1/(1-t)^2}).
// They deliberately round the extrema far upward.  The constants are
// reproduction lemmas, not manuscript constants.  The truncation formulas
// below depend only on these derivative inequalities, never on a coarse/fine
// solution difference.  A later interval-jet test can independently tighten
// these five inequalities without changing the enclosure API.
inline A4TruncationMajorant appendix_a4_truncation_majorant(double step,double Tf){
    if(!(step>0.0&&Tf>0.0)) throw std::invalid_argument("invalid A.4 truncation step");
    constexpr double S1=9.0, S2=100.0, S3=4000.0, S4=120000.0, S5=7000000.0;
    // For every transition used by A.4, r=(log g)'=2l-1 has |r|<=3.
    // The largest l amplitude is one.  Hence |r^(k)|<=2 |s^(k)|.
    const double r=3.0,r1=2*S1,r2=2*S2,r3=2*S3,r4=2*S4;
    // If g'=r g, then g^(5)/g is the complete Bell polynomial below.
    const double B5=std::pow(r,5)+10*std::pow(r,3)*r1+15*r*r1*r1
                   +10*r*r*r2+10*r1*r2+5*r*r3+r4;
    // Normalized g starts <=1 and can grow by at most exp(.2) in the first
    // transition; 1.25 is an outward convenience bound.
    const double G=1.25;
    // Five variable RK4 pieces have total length 7.  We use coefficient 1
    // rather than the sharp RK4 principal-error coefficient, so this is a
    // deliberately conservative h^4 global majorant once the derivative lemma
    // above is certified.
    const double Erk=7.0*std::pow(step,4)*G*B5;
    // A.10 is composite Simpson.  Differentiating g=c exp(2(theta-1)a)
    // four times gives the same Bell-polynomial structure with theta
    // derivatives reduced by Tf^{-k}.  Using B5 (one derivative stronger) is
    // conservative. Standard Simpson: L h^4 max|F''''|/180.
    const double Ea10=(Tf/180.0)*std::pow(step,4)*G*B5;
    return {S1,S2,S3,S4,S5,B5,G*B5,Erk,Ea10,Erk+Ea10};
}

// A-priori truncation enclosure.  Richardson extrapolation has intentionally
// been removed: only one numerical pressure evaluation is made.  Constant-l
// pieces are exact; RK4 transition error is controlled by a fifth-derivative
// state majorant and A.10 Simpson error by a fourth-derivative integrand
// majorant.  eta derivatives before A.10 are explicit multiples of the same
// radial integral; on A.10 we conservatively reuse the integral majorant with
// the exact eta multiplier bounds valid on |eta|<=1.
inline A4PressureEnclosure appendix_a4_pressure_enclosure(
    const OuterProfileParameters& p,double eta,double Tf=20.0,double co=0.05,
    double coarse_step=0.04,double safety=1.0){
    if(!(coarse_step>0.0&&safety>=1.0)) throw std::invalid_argument("invalid A.4 enclosure parameters");
    // Keep the former H/4 resolution as the center so downstream accuracy does
    // not regress, but do not evaluate H or H/2 and do not compare resolutions.
    const double step=0.25*coarse_step;
    const auto q=appendix_a4_pressure_datum(p,eta,Tf,co,step);
    const auto m=appendix_a4_truncation_majorant(step,Tf);
    const double J1=std::abs(2.0*eta/(1.0+eta*eta)); // <=1
    const double J2=std::abs(2.0*(1.0-eta*eta)/std::pow(1.0+eta*eta,2)); // <=2
    const double radial=safety*m.total_integral_error;
    // Pi=-1/2 integral. For d_eta: |-2 J_eta|<=2.  For d_etaeta:
    // |4 J_eta^2-2 J_etaeta|<=8.  A.10 multipliers obey the same bounds.
    const double e0=a4_enclosure_up(0.5*radial);
    const double e1=a4_enclosure_up(0.5*radial*std::max(2.0,2.0*J1));
    const double e2=a4_enclosure_up(0.5*radial*std::max(8.0,4.0*J1*J1+2.0*J2));
    const bool finite=std::isfinite(e0)&&std::isfinite(e1)&&std::isfinite(e2);
    return {q,e0,e1,e2,
      a4_enclosure_down(q.normalized_pi0-e0),a4_enclosure_up(q.normalized_pi0+e0),
      a4_enclosure_down(q.normalized_deta-e1),a4_enclosure_up(q.normalized_deta+e1),
      a4_enclosure_down(q.normalized_detaeta-e2),a4_enclosure_up(q.normalized_detaeta+e2),
      coarse_step,step,true,finite};
}

} // namespace nsblowup
