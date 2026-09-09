#pragma once

#include <cmath>
#include <stdexcept>

namespace nsblowup {

// Composite Simpson integration. The panel count is rounded up to an even value.
template <class Func>
double simpson_integral(Func&& f, double a, double b, int panels = 256) {
    if (panels < 2) {
        panels = 2;
    }
    if (panels % 2 != 0) {
        ++panels;
    }
    if (a == b) {
        return 0.0;
    }

    const double dx = (b - a) / static_cast<double>(panels);
    double sum = f(a) + f(b);
    for (int i = 1; i < panels; ++i) {
        const double x = a + dx * static_cast<double>(i);
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return sum * dx / 3.0;
}

// Equation (4.6): A_X(f) = X^{-1} integral_0^X f(x,eta) dx,
// with its smooth axis value A_X(f)(0,eta) = f(0,eta).
template <class Func>
double radial_average(Func&& f, double X, double eta, int panels = 256) {
    if (X < 0.0) {
        throw std::invalid_argument("X must be non-negative");
    }
    if (X == 0.0) {
        return f(0.0, eta);
    }
    return simpson_integral(
        [&](double x) { return f(x, eta); }, 0.0, X, panels) / X;
}

// First identity in equation (4.7). Since differentiation in eta commutes
// with the radial integral for a smooth profile, d_eta A_X(U) is computed as
// A_X(partial_eta U).
template <class UFunc, class DUDEtaFunc>
double leading_radial_flux_V0(
    UFunc&& U,
    DUDEtaFunc&& dU_deta,
    double X,
    double eta,
    double h,
    int panels = 256) {
    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }
    const double D = 0.5 - h;
    const double d = 1.0 - eta * eta;
    const double L = 1.0 - 2.0 * h * eta * eta;
    const double average_U = radial_average(U, X, eta, panels);
    const double d_eta_average_U = radial_average(dU_deta, X, eta, panels);

    return (X / L) * (
        2.0 * eta * U(X, eta)
        - 2.0 * D * eta * average_U
        - d * d_eta_average_U);
}

// For the regular-axis factorization E = sqrt(2X) F, equation (4.7) gives
// partial_X Pi = F^2. With Pi -> 0 at radial infinity (equation (4.25)),
// Pi(X,eta) = - integral_X^infinity F(x,eta)^2 dx.
// Numerically, infinity is replaced by X_max; callers must choose X_max in
// a region where the tail is negligible or exactly zero.
template <class FFunc>
double pressure_from_regular_F(
    FFunc&& F,
    double X,
    double eta,
    double X_max,
    int panels = 512) {
    if (X < 0.0 || X_max < X) {
        throw std::invalid_argument("require 0 <= X <= X_max");
    }
    return -simpson_integral(
        [&](double x) {
            const double value = F(x, eta);
            return value * value;
        },
        X,
        X_max,
        panels);
}

}  // namespace nsblowup
