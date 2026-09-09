#include "similarity_coordinates.hpp"
#include "similarity_operators.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

double profile(double X, double eta) {
    return std::exp(-X) * (1.0 + 0.2 * eta + 0.3 * eta * eta);
}

nsblowup::ProfileJet profile_jet(double X, double eta) {
    const double e = std::exp(-X);
    const double poly = 1.0 + 0.2 * eta + 0.3 * eta * eta;
    return {
        e * poly,
        e * (0.2 + 0.6 * eta),
        -X * e * poly,
    };
}

double physical_value(double r, double z, double tau, double h, double b) {
    const auto c = nsblowup::coordinates(r, z, tau, h);
    return std::pow(c.q, b) * profile(c.X, c.eta);
}

bool close(double a, double b, double rel = 2e-7) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel * scale;
}

int fail(const char* message, double numerical, double analytic) {
    std::cerr << message << ": numerical=" << numerical
              << " analytic=" << analytic << '\n';
    return 1;
}

}  // namespace

int main() {
    constexpr double h = 0.005;
    constexpr double D = 0.5 - h;
    constexpr double b = -0.37;
    constexpr double r = 0.41;
    constexpr double z = 0.07;
    constexpr double tau = 0.23;
    constexpr double eps = 1e-6;

    const auto c = nsblowup::coordinates(r, z, tau, h);
    const auto jet = profile_jet(c.X, c.eta);

    // t = 1 - tau, so increasing t decreases tau.
    const double dt_fd = (
        physical_value(r, z, tau - eps, h, b)
        - physical_value(r, z, tau + eps, h, b)
    ) / (2.0 * eps);
    const double dt_exact = std::pow(c.q, b - 1.0)
        * nsblowup::T_b(b, h, c.eta, jet);

    const double dz_fd = (
        physical_value(r, z + eps, tau, h, b)
        - physical_value(r, z - eps, tau, h, b)
    ) / (2.0 * eps);
    const double dz_exact = std::pow(c.q, b - D)
        * nsblowup::Z_b(b, h, c.eta, jet);

    if (!close(dt_fd, dt_exact)) {
        return fail("T_b finite-difference check failed", dt_fd, dt_exact);
    }
    if (!close(dz_fd, dz_exact)) {
        return fail("Z_b finite-difference check failed", dz_fd, dz_exact);
    }

    std::cout << "similarity-operator finite-difference checks passed\n";
    return 0;
}
