#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace nsblowup {

using BPolynomial = std::vector<double>;

inline BPolynomial appendix_b_Jnu(const BPolynomial& F, int nu) {
    if (nu != 1 && nu != 2) throw std::invalid_argument("Jnu requires nu=1 or 2");
    BPolynomial G(F.size() + 1, 0.0);
    for (std::size_t alpha = 0; alpha < F.size(); ++alpha) {
        const double a = static_cast<double>(alpha + 1);
        G[alpha + 1] = F[alpha] / (a * (static_cast<double>(alpha) + nu));
    }
    return G;
}

inline BPolynomial appendix_b_T(const BPolynomial& F, double chi, std::size_t n) {
    BPolynomial scaled = F;
    for (double& v : scaled) v *= 0.5 * chi;
    auto out = appendix_b_Jnu(scaled, 2);
    out.resize(n, 0.0);
    return out;
}

inline BPolynomial appendix_b_inverse_one_plus_T(
    const BPolynomial& rhs, double chi, std::size_t n) {
    BPolynomial out(n, 0.0);
    BPolynomial term = rhs;
    term.resize(n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        const double sign = (k % 2 == 0) ? 1.0 : -1.0;
        for (std::size_t j = 0; j < n; ++j) out[j] += sign * term[j];
        term = appendix_b_T(term, chi, n);
    }
    return out;
}

inline BPolynomial appendix_b_phi0_coefficients(double chi, std::size_t n) {
    BPolynomial c(n, 0.0);
    if (n == 0) return c;
    c[0] = 1.0;
    for (std::size_t a = 1; a < n; ++a) {
        c[a] = c[a - 1] * (-0.5 * chi) /
               (static_cast<double>(a) * static_cast<double>(a + 1));
    }
    return c;
}

} // namespace nsblowup
