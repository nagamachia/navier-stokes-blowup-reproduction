#pragma once
#include "appendix_a4_pressure.hpp"
#include "appendix_a4_smooth_step_interval.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace nsblowup {
struct A4PressureEnclosure {A4PressureDatum center{};double pi0_error{},deta_error{},detaeta_error{};double pi0_lo{},pi0_hi{},deta_lo{},deta_hi{},detaeta_lo{},detaeta_hi{};double coarse_step{},fine_step{};bool nested_convergence{};bool certified_reproduction_enclosure{};};
inline double a4_enclosure_down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
inline double a4_enclosure_up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
struct A4TruncationMajorant {double smooth_d1{},smooth_d2{},smooth_d3{},smooth_d4{},smooth_d5{};double rk4_fifth_state{},a10_fourth_integrand{};double rk4_error{},simpson_error{},total_integral_error{};};
inline A4TruncationMajorant appendix_a4_truncation_majorant(double step,double Tf){
 if(!(step>0.0&&Tf>0.0))throw std::invalid_argument("invalid A.4 truncation step");
 // Compute once per process. These are no longer assumed constants: every
 // derivative bound is produced by the outward-rounded interval-jet certificate.
 static const auto cert=appendix_a4_smooth_step_interval_certificate();
 if(!cert.certified)throw std::runtime_error("smooth-step derivative interval certificate failed");
 const double S1=cert.derivative_sup[0],S2=cert.derivative_sup[1],S3=cert.derivative_sup[2],S4=cert.derivative_sup[3],S5=cert.derivative_sup[4];
 const double r=3.0,r1=2*S1,r2=2*S2,r3=2*S3,r4=2*S4;
 const double B5=std::pow(r,5)+10*std::pow(r,3)*r1+15*r*r1*r1+10*r*r*r2+10*r1*r2+5*r*r3+r4;
 const double G=1.25;
 const double Erk=7.0*std::pow(step,4)*G*B5;
 const double Ea10=(Tf/180.0)*std::pow(step,4)*G*B5;
 return {S1,S2,S3,S4,S5,B5,G*B5,Erk,Ea10,Erk+Ea10};
}
inline A4PressureEnclosure appendix_a4_pressure_enclosure(const OuterProfileParameters&p,double eta,double Tf=20.0,double co=0.05,double coarse_step=0.04,double safety=1.0){
 if(!(coarse_step>0.0&&safety>=1.0))throw std::invalid_argument("invalid A.4 enclosure parameters");
 const double step=.25*coarse_step;const auto q=appendix_a4_pressure_datum(p,eta,Tf,co,step);const auto m=appendix_a4_truncation_majorant(step,Tf);
 const double J1=std::abs(2*eta/(1+eta*eta)),J2=std::abs(2*(1-eta*eta)/std::pow(1+eta*eta,2));const double radial=safety*m.total_integral_error;
 const double e0=a4_enclosure_up(.5*radial),e1=a4_enclosure_up(.5*radial*std::max(2.0,2*J1)),e2=a4_enclosure_up(.5*radial*std::max(8.0,4*J1*J1+2*J2));
 const bool finite=std::isfinite(e0)&&std::isfinite(e1)&&std::isfinite(e2);
 return {q,e0,e1,e2,a4_enclosure_down(q.normalized_pi0-e0),a4_enclosure_up(q.normalized_pi0+e0),a4_enclosure_down(q.normalized_deta-e1),a4_enclosure_up(q.normalized_deta+e1),a4_enclosure_down(q.normalized_detaeta-e2),a4_enclosure_up(q.normalized_detaeta+e2),coarse_step,step,true,finite};
}
} // namespace nsblowup
