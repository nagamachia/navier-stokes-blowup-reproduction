#include "appendix_b_parameter_scale.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : "appendix_b_diagnostics.csv";

    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;

    std::ofstream out(path);
    if (!out) {
        std::cerr << "failed to open output: " << path << '\n';
        return 1;
    }
    out << "U_tolerance,Ymax,max_log_abs_Zstar,eta_at_max,max_log_abs_u0,log_Lambda_min,log_Lambda_over_Pstar2\n";
    out << std::setprecision(17);
    for (double tol : {0.1, 0.05, 0.01}) {
        const auto a = nsblowup::appendix_b_lambda_scale_audit(
            p, tol, 4.1, 0.05, 20.0, 0.05, 0.04, 161);
        out << a.U_tolerance << ',' << a.Ymax << ','
            << a.max_log_abs_Zstar << ',' << a.eta_at_max << ','
            << a.max_log_abs_u0_Ymax << ','
            << a.log_Lambda_for_U_tolerance << ','
            << a.log_Lambda_over_Pstar2 << '\n';
    }
    return 0;
}
