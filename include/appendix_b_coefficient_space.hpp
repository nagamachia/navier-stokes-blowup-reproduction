#pragma once

#include <algorithm>
#include <cmath>
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

// Weighted analytic coefficient norm on |Y|<=rho.
inline double appendix_b_weighted_l1(const BPolynomial& F, double rho) {
    if (!(rho > 0.0)) throw std::invalid_argument("weighted norm requires rho>0");
    double s = 0.0;
    double w = 1.0;
    for (double c : F) {
        s += std::abs(c) * w;
        w *= rho;
    }
    return s;
}

// Exact coefficientwise operator bound in the weighted l1 norm:
// ||J_nu F||_rho <= rho/nu ||F||_rho, attained by the constant coefficient.
inline double appendix_b_Jnu_bound(double rho, int nu) {
    if (!(rho > 0.0) || (nu != 1 && nu != 2))
        throw std::invalid_argument("invalid Jnu bound parameters");
    return rho / static_cast<double>(nu);
}

inline double appendix_b_T_bound(double rho, double chi_abs_bound = 1.0) {
    if (!(chi_abs_bound >= 0.0)) throw std::invalid_argument("chi bound must be nonnegative");
    return 0.5 * chi_abs_bound * appendix_b_Jnu_bound(rho, 2);
}

struct AppendixBLinearOperatorAudit {
    double rho{};
    double J1_bound{};
    double J2_bound{};
    double T_bound{};
    double phi0_norm{};
    double inverse_image_norm{};
    double inverse_residual_norm{};
};

inline AppendixBLinearOperatorAudit appendix_b_linear_operator_audit(
    double chi, double rho = 4.1, std::size_t n = 32) {
    if (!(chi >= 0.0 && chi <= 1.0))
        throw std::invalid_argument("Appendix B chi must lie in [0,1]");
    BPolynomial one(n, 0.0);
    one[0] = 1.0;
    const auto inv = appendix_b_inverse_one_plus_T(one, chi, n);
    const auto Tinverse = appendix_b_T(inv, chi, n);
    BPolynomial residual(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        residual[i] = inv[i] + Tinverse[i] - one[i];

    AppendixBLinearOperatorAudit out;
    out.rho = rho;
    out.J1_bound = appendix_b_Jnu_bound(rho, 1);
    out.J2_bound = appendix_b_Jnu_bound(rho, 2);
    out.T_bound = appendix_b_T_bound(rho, chi);
    out.phi0_norm = appendix_b_weighted_l1(appendix_b_phi0_coefficients(chi, n), rho);
    out.inverse_image_norm = appendix_b_weighted_l1(inv, rho);
    out.inverse_residual_norm = appendix_b_weighted_l1(residual, rho);
    return out;
}

} // namespace nsblowup
