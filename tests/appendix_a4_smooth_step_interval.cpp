#include "appendix_a4_smooth_step_interval.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
int main(){
 const auto e=nsblowup::appendix_a4_smooth_step_endpoint_flatness();
 const auto c=nsblowup::appendix_a4_smooth_step_interval_certificate();
 std::cout<<std::setprecision(17)<<"endpoint_eps="<<e.eps<<" x_sup="<<e.x_sup<<" endpoint_certified="<<e.certified<<"\n";
 for(int k=0;k<5;++k)std::cout<<"endpoint_d"<<k+1<<"_sup="<<e.derivative_sup[k]<<" qd"<<k+1<<"_sup="<<e.q_derivative_sup[k]<<"\n";
 std::cout<<"boxes="<<c.boxes<<" endpoint_cut="<<c.endpoint_cut<<" certified="<<c.certified<<"\n";
 for(int k=0;k<5;++k)std::cout<<"d"<<k+1<<"_sup="<<c.derivative_sup[k]<<" claimed="<<c.claimed[k]<<"\n";
 if(!e.monotonicity_condition||!e.bell_polynomial_bound||!e.certified)return 1;
 if(!c.outward_rounded||!c.endpoint_flatness_certified||!c.certified)return 2;
 for(int k=0;k<5;++k)if(!(e.derivative_sup[k]>=0.0&&std::isfinite(e.derivative_sup[k])))return 3;
 return 0;
}
