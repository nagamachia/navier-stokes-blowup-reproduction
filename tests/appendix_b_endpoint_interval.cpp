#include "appendix_b_axis.hpp"
#include "appendix_b_endpoint_interval.hpp"

#include <cmath>
#include <iostream>

int main(){
    using namespace nsblowup;
    OuterProfileParameters p; p.lambda=1e-5; p.log_h=-80.0; p.log_Pstar=70.0;
    const auto sep=appendix_b_separation_audit(p,0.05,20.0,0.05,0.05,121);
    const auto q=appendix_b_endpoint_interval_audit(std::exp(p.log_h),0.05,sep.sigma_star,sep.eta0);
    if(!q.covers_eta_interval || !q.azimuthal_chi_certified || !q.axial_localized) return 1;
    if(!(q.axial_eta_lo<=sep.eta0 && sep.eta0<=q.axial_eta_hi)) return 2;
    if(!(q.min_abs_H_on_azimuthal_boxes>q.H_threshold)) return 3;
    std::cout<<"B.19 interval geometry: axial=["<<q.axial_eta_lo<<','<<q.axial_eta_hi
             <<"] width="<<q.axial_width<<" Hthr="<<q.H_threshold
             <<" boxes="<<q.azimuthal_boxes<<'/'<<q.axial_boxes<<'\n';
    return 0;
}
