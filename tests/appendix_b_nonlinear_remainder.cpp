#include "appendix_b_nonlinear_remainder.hpp"

#include <cmath>
#include <iostream>

int main() {
    using namespace nsblowup;
    constexpr std::size_t n = 8;

    AppendixBNonlinearRemainderInput in;
    in.eta = 0.2;
    in.h = 1e-6;
    in.Lambda = 1000.0;
    in.sigma_star = 0.01;
    in.g = 0.7;
    in.g_eta = -0.2;
    in.Phi = BPolynomial(n, 0.0);
    in.Phi_eta = BPolynomial(n, 0.0);
    in.u = BPolynomial(n, 0.0);
    in.u_eta = BPolynomial(n, 0.0);
    in.Phi[0] = 1.0;
    in.Phi[1] = -0.1;
    in.Phi_eta[0] = 0.03;
    in.u[1] = 0.4;
    in.u[2] = -0.05;
    in.u_eta[1] = -0.2;

    const auto out = appendix_b_nonlinear_remainder(in, n);
    for (const auto* p : {&out.B, &out.W, &out.Hc, &out.p, &out.p_eta, &out.R1, &out.R2}) {
        if (p->size() != n) return 1;
        for (double x : *p) if (!std::isfinite(x)) return 2;
    }

    // p=I(g^2 Phi^2): constant integrand g^2 gives Y coefficient g^2.
    if (std::abs(out.p[1] - in.g * in.g) > 1e-14) {
        std::cerr << "B.15 pressure primitive coefficient mismatch\n";
        return 3;
    }

    // p_eta=I(2 g g_eta Phi^2 + 2 g^2 Phi Phi_eta).
    const double expected_peta1 = 2.0 * in.g * in.g_eta +
                                  2.0 * in.g * in.g * in.Phi_eta[0];
    if (std::abs(out.p_eta[1] - expected_peta1) > 1e-14) {
        std::cerr << "B.15 pressure eta derivative mismatch\n";
        return 4;
    }

    // A_X preserves radial degree and divides coefficient alpha by alpha+1.
    const auto ax = appendix_b_AX(in.u, n);
    if (std::abs(ax[1] - in.u[1] / 2.0) > 1e-15 ||
        std::abs(ax[2] - in.u[2] / 3.0) > 1e-15) return 5;

    // D_X multiplies coefficient alpha by alpha.
    const auto dx = appendix_b_DX(in.u, n);
    if (std::abs(dx[1] - in.u[1]) > 1e-15 ||
        std::abs(dx[2] - 2.0 * in.u[2]) > 1e-15) return 6;

    return 0;
}
