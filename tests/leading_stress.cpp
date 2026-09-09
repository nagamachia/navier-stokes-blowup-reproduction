#include "leading_stress.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double rel = 5e-6) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel * scale;
}

int fail(const char* message, double numerical, double expected) {
    std::cerr << message << ": numerical=" << numerical
              << " expected=" << expected << '\n';
    return 1;
}

struct Manufactured {
    static double F(double X, double) { return std::exp(-X); }
    static double DXF(double X, double) { return -X * std::exp(-X); }
    static double dEtaF(double, double) { return 0.0; }

    static double U(double, double) { return 0.0; }
    static double DXU(double, double) { return 0.0; }
    static double dEtaU(double, double) { return 0.0; }

    static double Pi(double X, double) { return -0.5 * std::exp(-2.0 * X); }
    static double DXPi(double X, double) { return X * std::exp(-2.0 * X); }
    static double dEtaPi(double, double) { return 0.0; }
};

nsblowup::LeadingStressPoint point(double X, double eta, double h) {
    return nsblowup::leading_stress_point(
        Manufactured::F,
        Manufactured::DXF,
        Manufactured::dEtaF,
        Manufactured::U,
        Manufactured::DXU,
        Manufactured::dEtaU,
        Manufactured::Pi,
        Manufactured::DXPi,
        Manufactured::dEtaPi,
        X,
        eta,
        h,
        256);
}

}  // namespace

int main() {
    constexpr double h = 0.005;
    constexpr double A = 0.5 + h;

    // For F=e^{-X}, U=0, Pi=-e^{-2X}/2:
    // ell=1-X, W=1, Hc=D eta,
    // S_q=X-1-h,
    // S_n=2 eta e^{-2X}(X-A).
    for (double X : {0.0, 0.2, 0.7, 1.4}) {
        for (double eta : {-0.7, 0.0, 0.4}) {
            const auto p = point(X, eta, h);
            const double exact_Sq = X - 1.0 - h;
            const double exact_Sn = 2.0 * eta * std::exp(-2.0 * X) * (X - A);
            if (!close(p.ell, 1.0 - X, 1e-10)) {
                return fail("ell check failed", p.ell, 1.0 - X);
            }
            if (!close(p.S_q, exact_Sq, 1e-9)) {
                return fail("S_q check failed", p.S_q, exact_Sq);
            }
            if (!close(p.S_n, exact_Sn, 1e-9)) {
                return fail("S_n check failed", p.S_n, exact_Sn);
            }
        }
    }

    // Equation (4.10) smooth axis limits.
    const double eta = 0.4;
    const auto axis = point(0.0, eta, h);
    const double axis_Sq = -1.0 - h;
    const double axis_Sn = -2.0 * A * eta;
    if (!close(axis.Q_s, 0.5 * axis_Sq, 1e-8)) {
        return fail("Q_s axis limit failed", axis.Q_s, 0.5 * axis_Sq);
    }
    if (!close(axis.N_s, axis_Sn, 1e-8)) {
        return fail("N_s axis limit failed", axis.N_s, axis_Sn);
    }
    if (!close(axis.T0_theta, 0.0, 1e-8)
        || !close(axis.T0_z, 0.0, 1e-8)) {
        return fail("T0 axis regularity failed",
                    std::hypot(axis.T0_theta, axis.T0_z), 0.0);
    }

    // Numerically differentiate Q_s and N_s to verify the ODEs in (4.9):
    // D_X Q_s + (1+ell)Q_s = S_q,
    // D_X N_s + N_s = S_n.
    constexpr double X = 0.8;
    constexpr double eps = 2e-5;
    const auto pm = point(X - eps, eta, h);
    const auto p0 = point(X, eta, h);
    const auto pp = point(X + eps, eta, h);
    const double DX_Qs = X * (pp.Q_s - pm.Q_s) / (2.0 * eps);
    const double DX_Ns = X * (pp.N_s - pm.N_s) / (2.0 * eps);
    const double lhs_q = DX_Qs + (1.0 + p0.ell) * p0.Q_s;
    const double lhs_n = DX_Ns + p0.N_s;

    if (!close(lhs_q, p0.S_q, 2e-5)) {
        return fail("Q_s ODE identity failed", lhs_q, p0.S_q);
    }
    if (!close(lhs_n, p0.S_n, 2e-5)) {
        return fail("N_s ODE identity failed", lhs_n, p0.S_n);
    }

    std::cout << "leading-stress regression checks passed\n";
    return 0;
}
