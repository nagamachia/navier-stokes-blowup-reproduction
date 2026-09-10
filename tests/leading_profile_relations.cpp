#include "leading_profile_relations.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double rel = 2e-10) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel * scale;
}

int fail(const char* message, double numerical, double exact) {
    std::cerr << message << ": numerical=" << numerical
              << " exact=" << exact << '\n';
    return 1;
}

}  // namespace

int main() {
    constexpr double h = 0.005;
    constexpr double D = 0.5 - h;

    auto U = [](double X, double eta) {
        return (1.0 + eta * eta) * (1.0 + 2.0 * X);
    };
    auto dU_deta = [](double X, double eta) {
        return 2.0 * eta * (1.0 + 2.0 * X);
    };

    for (double X : {0.0, 0.1, 0.7, 1.3}) {
        for (double eta : {-0.8, -0.2, 0.0, 0.5}) {
            const double average = nsblowup::radial_average(U, X, eta);
            const double exact_average = (1.0 + eta * eta) * (1.0 + X);
            if (!close(average, exact_average)) {
                return fail("radial average check failed", average, exact_average);
            }

            const double d = 1.0 - eta * eta;
            const double L = 1.0 - 2.0 * h * eta * eta;
            const double exact_d_eta_average = 2.0 * eta * (1.0 + X);
            const double exact_V0 = (X / L) * (
                2.0 * eta * U(X, eta)
                - 2.0 * D * eta * exact_average
                - d * exact_d_eta_average);
            const double numerical_V0 = nsblowup::leading_radial_flux_V0(
                U, dU_deta, X, eta, h);
            if (!close(numerical_V0, exact_V0)) {
                return fail("V0 relation check failed", numerical_V0, exact_V0);
            }
        }
    }

    // E = sqrt(2X) exp(-X), hence F = exp(-X),
    // partial_X Pi = exp(-2X). For a finite upper limit X_max,
    // Pi = -(exp(-2X)-exp(-2X_max))/2 exactly.
    auto F = [](double X, double /*eta*/) { return std::exp(-X); };
    constexpr double X_max = 12.0;
    for (double X : {0.0, 0.2, 1.0, 3.0}) {
        const double numerical = nsblowup::pressure_from_regular_F(
            F, X, 0.3, X_max, 2048);
        const double exact = -0.5 * (
            std::exp(-2.0 * X) - std::exp(-2.0 * X_max));
        if (!close(numerical, exact, 2e-9)) {
            return fail("pressure integration check failed", numerical, exact);
        }
    }

    // Equation (4.15) manufactured polynomial check in regular variables.
    // U=a(1+X), F=b are chosen so all five moments are elementary.
    constexpr double eta = 0.25;
    constexpr double X = 1.4;
    const double a = 1.0 + eta * eta;
    const double b = 2.0 - eta;
    auto Um = [a](double x, double /*unused*/) { return a * (1.0 + x); };
    auto Fm = [b](double /*x*/, double /*unused*/) { return b; };
    const auto m = nsblowup::cumulative_radial_moments(Um, Fm, X, eta, 2048);

    const double exact_M = a * (X + X * X / 2.0);
    const double exact_I = b * X * X;
    const double exact_J = 2.0 * a * b * (X * X / 2.0 + X * X * X / 3.0);
    const double exact_S = a * a * (X + X * X + X * X * X / 3.0) - b * b * X * X / 2.0;
    const double exact_Cp = b * b * X;
    if (!close(m.M, exact_M, 2e-10)) return fail("M moment check failed", m.M, exact_M);
    if (!close(m.I, exact_I, 2e-10)) return fail("I moment check failed", m.I, exact_I);
    if (!close(m.J, exact_J, 2e-10)) return fail("J moment check failed", m.J, exact_J);
    if (!close(m.S, exact_S, 2e-10)) return fail("S moment check failed", m.S, exact_S);
    if (!close(m.Cp, exact_Cp, 2e-10)) return fail("Cp moment check failed", m.Cp, exact_Cp);

    std::cout << "leading-profile relation checks passed\n";
    return 0;
}
