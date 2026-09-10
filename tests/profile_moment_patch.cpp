#include "profile_moment_patch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

namespace {

double bump(double x, double a, double b) {
    if (x <= a || x >= b) return 0.0;
    const double y = (x - a) / (b - a);
    return std::exp(-1.0 / (y * (1.0 - y)));
}

bool close(double a, double b, double rel = 2e-7) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel * scale;
}

}  // namespace

int main() {
    constexpr double uc = 0.35;
    constexpr double e0 = 1.1;
    // The Jacobian block weights are 1, x^(alpha+1/2) for U and
    // x^(1/2), x^alpha, x^(alpha-1) for E. Avoid exponents that make
    // these powers coincide; alpha=0.1 is safely nondegenerate.
    constexpr double alpha = 0.1;

    nsblowup::ScalarProfile1d U0 = [](double) { return uc; };
    nsblowup::ScalarProfile1d E0 = [](double x) { return e0 * std::pow(x, alpha); };

    std::array<nsblowup::ScalarProfile1d, 2> ub{
        [](double x) { return bump(x, 1.1, 1.5); },
        [](double x) { return bump(x, 1.8, 2.3); }};
    std::array<nsblowup::ScalarProfile1d, 3> eb{
        [](double x) { return bump(x, 2.7, 3.2); },
        [](double x) { return bump(x, 3.7, 4.3); },
        [](double x) { return bump(x, 4.8, 5.6); }};

    nsblowup::FiveBumpMomentPatch patch(1.0, 5.8, U0, E0, ub, eb, 2048);
    const auto B = patch.jacobian_at_zero();

    // Test numerical invertibility by recovering a known coefficient vector,
    // rather than thresholding the raw determinant, whose magnitude depends
    // strongly on bump normalization and moment units.
    const nsblowup::Vector linear_target{0.2, -0.15, 0.1, -0.08, 0.06};
    nsblowup::Vector linear_rhs(5, 0.0);
    for (std::size_t i = 0; i < 5; ++i)
        for (std::size_t j = 0; j < 5; ++j) linear_rhs[i] += B[i][j] * linear_target[j];

    // Row-normalize before the solve, matching the nonlinear routine.
    nsblowup::Matrix scaled_B = B;
    nsblowup::Vector scaled_rhs = linear_rhs;
    for (std::size_t i = 0; i < 5; ++i) {
        double scale = 0.0;
        for (double value : B[i]) scale = std::max(scale, std::abs(value));
        if (!(scale > 0.0)) {
            std::cerr << "zero row in Corollary A.3 Jacobian\n";
            return 1;
        }
        for (double& value : scaled_B[i]) value /= scale;
        scaled_rhs[i] /= scale;
    }
    const auto linear_solved = nsblowup::solve_linear(scaled_B, scaled_rhs, 1e-14);
    for (std::size_t i = 0; i < linear_target.size(); ++i) {
        if (!close(linear_solved[i], linear_target[i], 5e-7)) {
            std::cerr << "Corollary A.3 linear Jacobian recovery failed at " << i << '\n';
            return 1;
        }
    }

    const nsblowup::Vector target{0.012, -0.009, 0.007, -0.006, 0.005};
    const auto discrepancy = patch.change(target).vector();
    const auto solved = patch.solve_for_change(discrepancy, 80);

    for (std::size_t i = 0; i < target.size(); ++i) {
        if (!close(solved[i], target[i])) {
            std::cerr << "five-bump coefficient recovery failed at " << i
                      << ": solved=" << solved[i] << " target=" << target[i] << '\n';
            return 1;
        }
    }

    const auto recovered = patch.change(solved).vector();
    for (std::size_t i = 0; i < recovered.size(); ++i) {
        if (!close(recovered[i], discrepancy[i], 5e-8)) {
            std::cerr << "five-moment discrepancy recovery failed at " << i << '\n';
            return 1;
        }
    }

    std::cout << "Corollary A.3 five-bump moment correction checks passed\n";
    return 0;
}
