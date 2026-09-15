#include "appendix_b_axis.hpp"
#include "appendix_b_endpoint_interval.hpp"
#include "appendix_b_Z_interval.hpp"
#include <cmath>
#include <iostream>
int main(){
 using namespace nsblowup;
 OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
 const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.02,241);
 const auto ep=appendix_b_endpoint_interval_audit(std::exp(p.log_h),0.05,sep.sigma_star,sep.eta0);
 if(!ep.azimuthal_chi_certified||!ep.axial_localized)return 1;
 const auto z=appendix_b_Z_interval_audit(p,ep.axial_eta_lo,ep.axial_eta_hi,sep.normalized_delta_star,0.05,20.0,0.05,0.01);
 std::cout<<"Appendix B axial interval: ["<<ep.axial_eta_lo<<','<<ep.axial_eta_hi<<"] center="<<z.Z_center
          <<" M="<<z.max_abs_Zprime<<" Z=["<<z.Z_lower<<','<<z.Z_upper<<"] minabs="<<z.min_abs_Z
          <<" delta="<<sep.normalized_delta_star<<" pressure_error="<<z.pressure_center_error
          <<" flags="<<z.pressure_derivatives_available<<'/'<<z.pressure_enclosure_certified<<'/'<<z.positive_on_interval<<'/'<<z.separated_from_delta<<'\n';
 if(!z.pressure_derivatives_available||!z.pressure_enclosure_certified||!z.positive_on_interval||!z.separated_from_delta)return 2;
 if(!(z.min_abs_Z>sep.normalized_delta_star)||!(z.pressure_center_error>=0.0))return 3;
 return 0;
}
