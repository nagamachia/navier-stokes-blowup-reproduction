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
    constexpr double alpha = 0.5;  // one of the admissible Corollary A.3 exponents

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
    const double det = nsblowup::determinant(B, 1e-18);
    if (!(std::abs(det) > 1e-13)) {
        std::cerr << "Corollary A.3 Jacobian is numerically singular: det=" << det << '\n';
        return 1;
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
