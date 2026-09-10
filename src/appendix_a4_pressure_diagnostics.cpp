#include "appendix_a4_pressure.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : "appendix_a4_pressure.csv";
    std::ofstream f(out);
    if (!f) {
        std::cerr << "cannot open output: " << out << '\n';
        return 1;
    }

    nsblowup::OuterProfileParameters p;
    p.Md = 4.0;
    p.log_Pstar = 70.0;
    p.lambda = 1.0e-5;
    p.log_h = -80.0;
    p.validate();

    constexpr double Tf = 20.0;
    constexpr double co = 0.05;
    constexpr double step = 0.03;

    f << "eta,lambda,log_Pstar,log_h,pi0_over_Pstar2,dpi0_deta_over_Pstar2,"
         "d2pi0_deta2_over_Pstar2,inner_bound_over_Pstar2,a22_margin,"
         "eta_times_dpi0_deta_over_Pstar2,pi0,pressure_y0_over_Pstar2,"
         "exterior_hold_length,terminal_Qp\n";
    f << std::setprecision(17);

    for (double eta : {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto d = nsblowup::appendix_a4_pressure_datum(p, eta, Tf, co, step);
        const double pi0 = nsblowup::appendix_a4_scale_pressure(d.normalized_pi0, p.log_Pstar);
        const double p0 = nsblowup::appendix_a4_inner_pressure_normalized(d, 0.0);
        f << eta << ',' << p.lambda << ',' << p.log_Pstar << ',' << p.log_h << ','
          << d.normalized_pi0 << ',' << d.normalized_deta << ','
          << d.normalized_detaeta << ',' << d.inner_contribution << ','
          << (d.inner_contribution - d.normalized_pi0) << ','
          << eta*d.normalized_deta << ',' << pi0 << ',' << p0 << ','
          << d.exterior_hold_length << ',' << d.terminal_Qp << '\n';
    }

    std::cout << "wrote Appendix A.4 pressure diagnostics to " << out << '\n';
    return 0;
}
