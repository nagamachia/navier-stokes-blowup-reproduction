#include "appendix_a3_global.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : "appendix_a3_closure.csv";
    std::ofstream f(out);
    if (!f) {
        std::cerr << "cannot open output: " << out << '\n';
        return 1;
    }

    f << "lambda,eta,Kb,amplitude,m_pulse,j_pulse,s_pulse,pre_M_bump,pre_J_bump,"
         "c1_MJ,c2_MJ,residual_M,residual_J,c1_IP,c2_IP,residual_I,residual_pressure,Qp,Qs0\n";
    f << std::setprecision(17);

    const double Kb = nsblowup::appendix_a3_Kb();
    for (double lambda : {0.08, 0.05, 0.02}) {
        nsblowup::OuterProfileParameters p;
        p.Md = 2.0;
        p.log_Pstar = 25.0;
        p.lambda = lambda;
        p.log_h = -40.0;
        p.validate();

        for (double eta : {-0.5, 0.0, 0.5}) {
            const double amp = nsblowup::appendix_a3_solve_amplitude(eta);
            const auto pre = nsblowup::appendix_a3_integrate_to_pulse(p, eta, 0.03);
            const auto mj = nsblowup::appendix_a3_close_MJ_from_schedule(p, eta, amp, 0.03);

            // The angular correction is still evaluated as its exact local
            // A.11 quadratic map.  Wiring its discrepancy to the full
            // post-pulse schedule is the next integration step.
            const double dI = std::pow(lambda, 28.0);
            const auto ip = nsblowup::appendix_a3_close_angular_I_pressure(lambda, dI, 0.0);

            const double h = 1.0e-3;
            const double rho_o = 1.0e-4;
            const double Qp = nsblowup::outer_terminal_Qp(h, rho_o, 800);
            const double Qs0 = nsblowup::appendix_a3_terminal_Qs(0.0, h, rho_o, 800);

            f << lambda << ',' << eta << ',' << Kb << ',' << amp << ','
              << pre.m_at_pulse << ',' << pre.j_at_pulse << ',' << pre.s_at_pulse << ','
              << pre.pre_M_at_first_bump << ',' << pre.pre_J_at_first_bump << ','
              << mj.c1 << ',' << mj.c2 << ',' << mj.residual_M << ',' << mj.residual_J << ','
              << ip.c1 << ',' << ip.c2 << ',' << ip.residual_I << ',' << ip.residual_pressure << ','
              << Qp << ',' << Qs0 << '\n';
        }
    }

    std::cout << "wrote Appendix A.3 diagnostics to " << out << '\n';
    return 0;
}
