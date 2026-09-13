#include "appendix_b_axis.hpp"
#include "appendix_b_coefficient_space.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {
double eval(const nsblowup::BPolynomial& c, double y) {
    double s = 0.0;
    for (auto it = c.rbegin(); it != c.rend(); ++it) s = s * y + *it;
    return s;
}
}

int main() {
    nsblowup::OuterProfileParameters p;
    p.lambda = 1e-5;
    p.log_h = -80.0;
    p.log_Pstar = 70.0;
    const auto sep = nsblowup::appendix_b_separation_audit(p, 0.05, 20.0, 0.05, 0.04, 161);
    const double h = std::exp(p.log_h);

    for (double eta : {-1.0, -0.5, sep.eta0, 0.0, 0.5, 1.0}) {
        const double H = nsblowup::appendix_b_Hstar(eta, h, 0.05);
        const double chi = H * H / (H * H + sep.sigma_star * sep.sigma_star);
        constexpr std::size_t n = 24;
        nsblowup::BPolynomial one(n, 0.0);
        one[0] = 1.0;
        const auto inv = nsblowup::appendix_b_inverse_one_plus_T(one, chi, n);
        const auto exact = nsblowup::appendix_b_phi0_coefficients(chi, n);

        double coeff_err = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            coeff_err = std::max(coeff_err, std::abs(inv[i] - exact[i]));
        if (!(coeff_err < 1e-13)) {
            std::cerr << "B.2 coefficient inverse mismatch\n";
            return 1;
        }

        const auto Tinverse = nsblowup::appendix_b_T(inv, chi, n);
        double residual = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            residual = std::max(residual, std::abs(inv[i] + Tinverse[i] - one[i]));
        if (!(residual < 1e-13)) {
            std::cerr << "B.2 (1+T) inverse residual too large\n";
            return 2;
        }

        for (int k = 0; k <= 41; ++k) {
            const double Y = 4.1 * static_cast<double>(k) / 41.0;
            if (!(std::abs(eval(inv, Y) - nsblowup::appendix_b_f0(Y * chi)) < 1e-10)) {
                std::cerr << "B.2 coefficient profile disagrees with f0\n";
                return 3;
            }
        }
    }

    const nsblowup::BPolynomial constant{2.0};
    const auto j1 = nsblowup::appendix_b_Jnu(constant, 1);
    const auto j2 = nsblowup::appendix_b_Jnu(constant, 2);
    if (!(std::abs(j1[1] - 2.0) < 1e-15 && std::abs(j2[1] - 1.0) < 1e-15)) {
        std::cerr << "B.5 Jnu coefficient formula failed\n";
        return 4;
    }

    return 0;
}
