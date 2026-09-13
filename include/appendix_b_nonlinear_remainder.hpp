#pragma once

#include "appendix_b_axis.hpp"
#include "appendix_b_coefficient_space.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace nsblowup {

inline BPolynomial appendix_b_poly_add(const BPolynomial& a, const BPolynomial& b, std::size_t n) {
    BPolynomial out(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (i < a.size()) out[i] += a[i];
        if (i < b.size()) out[i] += b[i];
    }
    return out;
}

inline BPolynomial appendix_b_poly_scale(const BPolynomial& a, double s, std::size_t n) {
    BPolynomial out(n, 0.0);
    for (std::size_t i = 0; i < std::min(n, a.size()); ++i) out[i] = s * a[i];
    return out;
}

inline BPolynomial appendix_b_poly_mul(const BPolynomial& a, const BPolynomial& b, std::size_t n) {
    BPolynomial out(n, 0.0);
    for (std::size_t i = 0; i < a.size() && i < n; ++i)
        for (std::size_t j = 0; j < b.size() && i + j < n; ++j)
            out[i + j] += a[i] * b[j];
    return out;
}

// D_X = Y d/dY after Y=Lambda X.
inline BPolynomial appendix_b_DX(const BPolynomial& a, std::size_t n) {
    BPolynomial out(n, 0.0);
    for (std::size_t i = 1; i < std::min(n, a.size()); ++i)
        out[i] = static_cast<double>(i) * a[i];
    return out;
}

// A_X(F)=Y^{-1} int_0^Y F(Y') dY'.
inline BPolynomial appendix_b_AX(const BPolynomial& a, std::size_t n) {
    BPolynomial out(n, 0.0);
    for (std::size_t i = 0; i < std::min(n, a.size()); ++i)
        out[i] = a[i] / static_cast<double>(i + 1);
    return out;
}

// I(F)=int_0^Y F(Y') dY'.
inline BPolynomial appendix_b_I(const BPolynomial& a, std::size_t n) {
    BPolynomial out(n, 0.0);
    for (std::size_t i = 0; i < a.size() && i + 1 < n; ++i)
        out[i + 1] = a[i] / static_cast<double>(i + 1);
    return out;
}

struct AppendixBNonlinearRemainderInput {
    double eta{};
    double h{};
    double j0{0.05};
    double Lambda{};
    double sigma_star{};
    double g{};       // g=phi*/C at this eta
    double g_eta{};   // eta derivative of g
    BPolynomial Phi;
    BPolynomial Phi_eta;
    BPolynomial u;
    BPolynomial u_eta;
};

struct AppendixBNonlinearRemainder {
    BPolynomial B;
    BPolynomial W;
    BPolynomial Hc;
    BPolynomial p;
    BPolynomial p_eta;
    BPolynomial R1;
    BPolynomial R2;
};

// Exact radial-coefficient transcription of the remainders following (B.15).
// eta differentiation is supplied by the caller; no finite-difference formula is hidden here.
// Products are exact modulo the requested radial truncation n.
inline AppendixBNonlinearRemainder appendix_b_nonlinear_remainder(
    const AppendixBNonlinearRemainderInput& in, std::size_t n) {
    if (!(in.eta >= -1.0 && in.eta <= 1.0) || !(in.h > 0.0) ||
        !(in.Lambda > 0.0) || !(in.sigma_star > 0.0) || n == 0)
        throw std::invalid_argument("invalid Appendix B nonlinear remainder input");

    const double eta = in.eta;
    const double h = in.h;
    const double A = appendix_b_A(h);
    const double D = appendix_b_D(h);
    const double d = appendix_b_d(eta);
    const double L = appendix_b_L(eta, h);
    const double Ustar = appendix_b_Ustar(eta, in.j0);
    const double Ustar_eta = 4.0;
    const double Hstar = appendix_b_Hstar(eta, h, in.j0);
    const double Wstar = 1.0 - d * Ustar_eta - 2.0 * D * eta * Ustar;
    const double zeta = -L * Hstar / (Hstar * Hstar + in.sigma_star * in.sigma_star);

    const auto AXu = appendix_b_AX(in.u, n);
    const auto AXu_eta = appendix_b_AX(in.u_eta, n);

    auto B = appendix_b_poly_add(
        appendix_b_poly_scale(AXu, -2.0 * D * eta, n),
        appendix_b_poly_scale(AXu_eta, -d, n), n);

    BPolynomial one(n, 0.0); one[0] = 1.0;
    auto W = appendix_b_poly_add(
        appendix_b_poly_scale(one, Wstar, n),
        appendix_b_poly_scale(B, 1.0 / in.Lambda, n), n);
    auto Hc = appendix_b_poly_add(
        appendix_b_poly_scale(one, Hstar, n),
        appendix_b_poly_scale(in.u, d / in.Lambda, n), n);

    const auto Phi2 = appendix_b_poly_mul(in.Phi, in.Phi, n);
    auto p = appendix_b_I(appendix_b_poly_scale(Phi2, in.g * in.g, n), n);
    auto p_eta_integrand = appendix_b_poly_add(
        appendix_b_poly_scale(Phi2, 2.0 * in.g * in.g_eta, n),
        appendix_b_poly_scale(appendix_b_poly_mul(in.Phi, in.Phi_eta, n),
                              2.0 * in.g * in.g, n), n);
    auto p_eta = appendix_b_I(p_eta_integrand, n);

    // U = U* + Lambda^{-1} u.
    auto U = appendix_b_poly_add(
        appendix_b_poly_scale(one, Ustar, n),
        appendix_b_poly_scale(in.u, 1.0 / in.Lambda, n), n);
    auto one_minus_2etaU = appendix_b_poly_add(
        one, appendix_b_poly_scale(U, -2.0 * eta, n), n);

    // L R1 = [W + h(1-2 eta U) + d u zeta*] Phi
    //        + W D_X Phi + Hc Phi_eta.
    auto bracket1 = appendix_b_poly_add(
        W, appendix_b_poly_scale(one_minus_2etaU, h, n), n);
    bracket1 = appendix_b_poly_add(
        bracket1, appendix_b_poly_scale(in.u, d * zeta, n), n);
    auto LR1 = appendix_b_poly_mul(bracket1, in.Phi, n);
    LR1 = appendix_b_poly_add(LR1,
        appendix_b_poly_mul(W, appendix_b_DX(in.Phi, n), n), n);
    LR1 = appendix_b_poly_add(LR1,
        appendix_b_poly_mul(Hc, in.Phi_eta, n), n);
    auto R1 = appendix_b_poly_scale(LR1, 1.0 / L, n);

    // L R2 = [A(1-4 eta U*) + d U*_eta]u - 2 A eta Lambda^{-1}u^2
    //        + W D_X u + H* u_eta + d Lambda^{-1}u u_eta
    //        - 4 A eta p + d p_eta - 2 eta D_X p.
    const double linear_u = A * (1.0 - 4.0 * eta * Ustar) + d * Ustar_eta;
    auto LR2 = appendix_b_poly_scale(in.u, linear_u, n);
    LR2 = appendix_b_poly_add(LR2,
        appendix_b_poly_scale(appendix_b_poly_mul(in.u, in.u, n),
                              -2.0 * A * eta / in.Lambda, n), n);
    LR2 = appendix_b_poly_add(LR2,
        appendix_b_poly_mul(W, appendix_b_DX(in.u, n), n), n);
    LR2 = appendix_b_poly_add(LR2,
        appendix_b_poly_scale(in.u_eta, Hstar, n), n);
    LR2 = appendix_b_poly_add(LR2,
        appendix_b_poly_scale(appendix_b_poly_mul(in.u, in.u_eta, n),
                              d / in.Lambda, n), n);
    LR2 = appendix_b_poly_add(LR2, appendix_b_poly_scale(p, -4.0 * A * eta, n), n);
    LR2 = appendix_b_poly_add(LR2, appendix_b_poly_scale(p_eta, d, n), n);
    LR2 = appendix_b_poly_add(LR2,
        appendix_b_poly_scale(appendix_b_DX(p, n), -2.0 * eta, n), n);
    auto R2 = appendix_b_poly_scale(LR2, 1.0 / L, n);

    return {B, W, Hc, p, p_eta, R1, R2};
}

} // namespace nsblowup
