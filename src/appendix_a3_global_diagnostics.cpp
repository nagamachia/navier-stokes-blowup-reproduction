#include "appendix_a3_global_closure.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : "appendix_a3_global.csv";
    std::ofstream f(out);
    if (!f) {
        std::cerr << "cannot open output: " << out << '\n';
        return 1;
    }

    nsblowup::OuterProfileParameters p;
    p.Md = 2.0;
    p.log_Pstar = 25.0;
    p.lambda = 1.0e-5;
    p.log_h = -40.0;
    p.validate();

    const double Kb = nsblowup::appendix_a3_Kb();
    const double principal_constant = (1.0 - std::exp(-26.0)) / 4.0;
    const double step = 0.05;

    f << "eta,lambda,amplitude,s_infinity,lambda_s_infinity,a19_principal,a19_remainder,"
         "m_pulse,j_pulse,s_pulse,rI_pulse,pulse_residual_M,pulse_residual_J,"
         "rI_exterior_error,angular_c1,angular_c2,exterior_hold_length,terminal_weight\n";
    f << std::setprecision(17);

    for (double eta : {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto pre = nsblowup::appendix_a3_integrate_to_pulse(p, eta, step);
        const auto g = nsblowup::appendix_a3_solve_global_amplitude_direct(
            p, eta, 20.0, 0.05, step, 34);
        const double principal = g.amplitude * g.amplitude * Kb - principal_constant;
        const double scaled = p.lambda * g.s_infinity;
        const double remainder = scaled - principal;

        f << eta << ',' << p.lambda << ',' << g.amplitude << ','
          << g.s_infinity << ',' << scaled << ',' << principal << ',' << remainder << ','
          << pre.m_at_pulse << ',' << pre.j_at_pulse << ',' << pre.s_at_pulse << ',' << pre.rI_at_pulse << ','
          << g.pulse_residual_M << ',' << g.pulse_residual_J << ','
          << (g.rI_exterior_start - g.rI_target) << ','
          << g.angular_c1 << ',' << g.angular_c2 << ','
          << g.exterior_hold_length << ',' << g.terminal_weight << '\n';
    }

    std::cout << "wrote global Appendix A.3 diagnostics to " << out << '\n';
    return 0;
}
