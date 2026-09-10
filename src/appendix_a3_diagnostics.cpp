#include "appendix_a3_closure.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : "appendix_a3_closure.csv";
    std::ofstream f(out);
    if (!f) {
        std::cerr << "cannot open output: " << out << '\n';
        return 1;
    }

    f << "lambda,Kb,amplitude,c1_MJ,c2_MJ,residual_M,residual_J,"
         "c1_IP,c2_IP,residual_I,residual_pressure,Qp,Qs0\n";
    f << std::setprecision(17);

    const double Kb = nsblowup::appendix_a3_Kb();
    for (double lambda : {0.08, 0.05, 0.02}) {
        const double amp = nsblowup::appendix_a3_solve_amplitude(0.0);
        // A.14/A.15 predict exponentially small normalized discrepancies at
        // the late pulse bumps. Use a representative exp(-1/lambda) scale.
        const double eps = std::exp(-1.0 / lambda);
        const auto mj = nsblowup::appendix_a3_close_pulse_MJ(
            lambda, eps, -0.7 * eps, amp);

        // A.3 gives the angular correction discrepancy O(lambda^28).
        const double dI = std::pow(lambda, 28.0);
        const auto ip = nsblowup::appendix_a3_close_angular_I_pressure(
            lambda, dI, 0.0);

        const double h = 1.0e-3;
        const double rho_o = 1.0e-4;
        const double Qp = nsblowup::outer_terminal_Qp(h, rho_o, 800);
        const double Qs0 = nsblowup::appendix_a3_terminal_Qs(0.0, h, rho_o, 800);

        f << lambda << ',' << Kb << ',' << amp << ','
          << mj.c1 << ',' << mj.c2 << ',' << mj.residual_M << ',' << mj.residual_J << ','
          << ip.c1 << ',' << ip.c2 << ',' << ip.residual_I << ',' << ip.residual_pressure << ','
          << Qp << ',' << Qs0 << '\n';
    }

    std::cout << "wrote Appendix A.3 diagnostics to " << out << '\n';
    return 0;
}
