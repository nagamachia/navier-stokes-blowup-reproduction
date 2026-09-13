#include "appendix_a7_heat_compensation.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double atol, double rtol) {
    return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}

}

int main() {
    using namespace nsblowup;

    const double lambda = 0.03;
    for (double eta : {-0.9, -0.5, 0.0, 0.5, 0.9}) {
        AppendixA7HeatCompensation patch(eta, lambda);
        const Matrix B = patch.jacobian_at_zero();

        // Proposition A.7 / Lemma A.1: the three distinct power weights
        // on three ordered supports give an invertible 3x3 Jacobian.
        Vector scale(3, 0.0);
        Matrix scaled = B;
        for (std::size_t i = 0; i < 3; ++i) {
            for (double v : B[i]) scale[i] = std::max(scale[i], std::abs(v));
            for (double& v : scaled[i]) v /= scale[i];
        }
        const double det = determinant(scaled);
        if (!(std::abs(det) > 1e-4)) {
            std::cerr << "A.7 normalized Jacobian nearly singular at eta=" << eta
                      << ": det=" << det << "\n";
            return 1;
        }

        const Vector known{1.2e-4, -8.0e-5, 5.0e-5};
        const Vector target = patch.change(known).vector();
        const Vector solved = patch.solve_for_change(target);
        const Vector recovered = patch.change(solved).vector();

        for (std::size_t i = 0; i < 3; ++i) {
            if (!close(recovered[i], target[i], 2e-12, 2e-10)) {
                std::cerr << "A.7 moment recovery failed at eta=" << eta
                          << " row=" << i << "\n";
                return 1;
            }
            if (!close(solved[i], known[i], 5e-10, 5e-6)) {
                std::cerr << "A.7 coefficient recovery failed at eta=" << eta
                          << " col=" << i << " solved=" << solved[i]
                          << " known=" << known[i] << "\n";
                return 1;
            }
        }

        // U is identically zero throughout this correction, so M and J are
        // structurally unchanged; the implemented map contains only Cp,S,I.
        if (target.size() != 3) {
            std::cerr << "unexpected A.7 correction dimension\n";
            return 1;
        }
    }

    // f(eta)=(1+eta^2)^-1 is even, so the finite-dimensional map must be even.
    AppendixA7HeatCompensation plus(0.63, lambda), minus(-0.63, lambda);
    const Vector c{7e-5, -3e-5, 4e-5};
    const Vector vp = plus.change(c).vector();
    const Vector vm = minus.change(c).vector();
    for (std::size_t i = 0; i < 3; ++i) {
        if (!close(vp[i], vm[i], 1e-12, 1e-12)) {
            std::cerr << "A.7 eta symmetry failed in row " << i << "\n";
            return 1;
        }
    }

    std::cout << "Appendix A.7 heat compensation regression passed\n";
    return 0;
}
