#include "leading_field.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double rel = 2e-12) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel * scale;
}

int fail(const char* message, double numerical, double expected) {
    std::cerr << message << ": numerical=" << numerical
              << " expected=" << expected << '\n';
    return 1;
}

}  // namespace

int main() {
    constexpr double h = 0.005;
    constexpr double A = 0.5 + h;
    constexpr double q = 0.37;
    constexpr double eta = 0.2;

    auto F = [](double X, double e) { return 1.0 + 0.1 * X + 0.2 * e; };
    auto U = [](double X, double e) { return 0.7 - 0.2 * X + 0.1 * e; };
    auto v0 = [](double X, double e) { return 2.0 + 0.3 * X - 0.1 * e; };
    auto Pi = [](double X, double e) { return -0.4 + 0.05 * X + 0.02 * e; };

    const auto axis = nsblowup::cartesian_leading_field(
        F, U, v0, Pi, 0.0, 0.0, q, eta, h);
    if (!close(axis.u1, 0.0) || !close(axis.u2, 0.0)) {
        return fail("transverse axis velocity must vanish",
                    std::hypot(axis.u1, axis.u2), 0.0);
    }
    if (!close(axis.u3, std::pow(q, -A) * U(0.0, eta))) {
        return fail("axial axis velocity mismatch",
                    axis.u3, std::pow(q, -A) * U(0.0, eta));
    }
    if (!close(axis.pressure, std::pow(q, -2.0 * A) * Pi(0.0, eta))) {
        return fail("axis pressure mismatch",
                    axis.pressure, std::pow(q, -2.0 * A) * Pi(0.0, eta));
    }

    for (double r : {1e-8, 1e-5, 1e-3, 0.1}) {
        constexpr double theta = 0.73;
        const double x1 = r * std::cos(theta);
        const double x2 = r * std::sin(theta);
        const double X = r * r / (2.0 * q);
        const auto field = nsblowup::cartesian_leading_field(
            F, U, v0, Pi, x1, x2, q, eta, h);

        const double u_r = (x1 * field.u1 + x2 * field.u2) / r;
        const double u_theta = (-x2 * field.u1 + x1 * field.u2) / r;
        const double exact_u_r = r * v0(X, eta) / (2.0 * q);
        const double exact_u_theta = std::pow(q, -A - 0.5) * r * F(X, eta);

        if (!close(u_r, exact_u_r, 5e-12)) {
            return fail("radial reconstruction failed", u_r, exact_u_r);
        }
        if (!close(u_theta, exact_u_theta, 5e-12)) {
            return fail("azimuthal reconstruction failed", u_theta, exact_u_theta);
        }
    }

    // Direction independence as r -> 0: transverse velocity is linear in x.
    constexpr double r = 1e-7;
    const auto f0 = nsblowup::cartesian_leading_field(
        F, U, v0, Pi, r, 0.0, q, eta, h);
    const auto f1 = nsblowup::cartesian_leading_field(
        F, U, v0, Pi, 0.0, r, q, eta, h);
    if (!(std::hypot(f0.u1, f0.u2) < 1e-5
          && std::hypot(f1.u1, f1.u2) < 1e-5)) {
        return fail("near-axis transverse velocity is not O(r)",
                    std::max(std::hypot(f0.u1, f0.u2), std::hypot(f1.u1, f1.u2)),
                    0.0);
    }

    std::cout << "leading-field axis regularity checks passed\n";
    return 0;
}
