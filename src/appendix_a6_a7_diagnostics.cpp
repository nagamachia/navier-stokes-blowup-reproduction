#include "appendix_a6_heat.hpp"
#include "appendix_a7_heat_compensation.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    using namespace nsblowup;
    const std::string out = argc > 1 ? argv[1] : "appendix_a6_a7_heat.csv";
    std::ofstream os(out);
    if (!os) {
        std::cerr << "failed to open " << out << "\n";
        return 1;
    }

    os << "eta,heat_ratio_x100,heat_ode_residual,jacobian_det_scaled,"
          "recovery_error_inf,c0,c1,c2\n";
    os << std::setprecision(15);

    constexpr double h = 0.08;
    constexpr double lambda = 0.03;
    const Vector known{1.2e-4, -8.0e-5, 5.0e-5};

    for (int k = -9; k <= 9; ++k) {
        const double eta = static_cast<double>(k) / 10.0;
        const double Z = 2.0 * (1.0 - eta * eta) / 100.0;
        const double heat_ratio = appendix_a6_heat_ratio(100.0, eta, h);
        const double ode = appendix_a6_heat_ode_residual(Z, h);

        AppendixA7HeatCompensation patch(eta, lambda);
        Matrix B = patch.jacobian_at_zero();
        for (std::size_t i = 0; i < B.size(); ++i) {
            double scale = 0.0;
            for (double v : B[i]) scale = std::max(scale, std::abs(v));
            for (double& v : B[i]) v /= scale;
        }
        const double det = determinant(B);
        const Vector target = patch.change(known).vector();
        const Vector solved = patch.solve_for_change(target);
        const Vector recovered = patch.change(solved).vector();
        double err = 0.0;
        for (std::size_t i = 0; i < target.size(); ++i)
            err = std::max(err, std::abs(recovered[i] - target[i]));

        os << eta << ',' << heat_ratio << ',' << ode << ',' << det << ','
           << err << ',' << solved[0] << ',' << solved[1] << ',' << solved[2] << '\n';
    }

    std::cout << "wrote Appendix A.6/A.7 diagnostics to " << out << "\n";
    return 0;
}
