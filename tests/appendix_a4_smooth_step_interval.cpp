#include "appendix_a4_smooth_step_interval.hpp"
#include <iostream>
#include <iomanip>
int main(){
 const auto c=nsblowup::appendix_a4_smooth_step_interval_certificate();
 std::cout<<std::setprecision(17)<<"boxes="<<c.boxes<<" endpoint_cut="<<c.endpoint_cut<<" certified="<<c.certified<<"\n";
 for(int k=0;k<5;++k)std::cout<<"d"<<k+1<<"_sup="<<c.derivative_sup[k]<<" claimed="<<c.claimed[k]<<"\n";
 if(!c.outward_rounded||!c.certified)return 1;
 return 0;
}
