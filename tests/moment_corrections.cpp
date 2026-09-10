#include "moment_corrections.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

double bump(double x, double a, double b) {
    if (x <= a || x >= b) return 0.0;
    const double y = (x - a) / (b - a);
    return std::exp(-1.0 / (y * (1.0 - y)));
}

bool close(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) <= tol * std::max({1.0, std::abs(a), std::abs(b)});
}

int fail(const char* message) {
    std::cerr << message << '\n';
    return 1;
}

}  // namespace

int main() {
    const nsblowup::Vector alphas{-1.0, 0.5, 2.0};
    const std::vector<std::pair<double, double>> intervals{{0.5, 0.9}, {1.2, 1.7}, {2.1, 2.8}};
    std::vector<std::function<double(double)>> bumps;
    for (const auto [a, b] : intervals) {
        bumps.push_back([a, b](double x) { return bump(x, a, b); });
    }

    const auto B = nsblowup::power_moment_matrix(alphas, intervals, bumps, 2048);
    const double det = nsblowup::determinant(B);
    if (!(std::abs(det) > 1e-12)) return fail("Lemma A.1 moment matrix should be invertible");

    const nsblowup::Vector exact{0.12, -0.08, 0.05};
    nsblowup::Vector rhs(B.size(), 0.0);
    for (std::size_t i = 0; i < B.size(); ++i)
        for (std::size_t j = 0; j < exact.size(); ++j) rhs[i] += B[i][j] * exact[j];
    const auto solved = nsblowup::solve_linear(B, rhs);
    for (std::size_t i = 0; i < exact.size(); ++i)
        if (!close(solved[i], exact[i], 2e-10)) return fail("linear moment correction recovery failed");

    // A small nonlinear system of the Lemma A.2 form F(c)=Bc+Q(c,c)=d.
    const nsblowup::Matrix B2{{1.2, 0.1}, {-0.05, 0.9}};
    const nsblowup::BilinearMap Q = [](const nsblowup::Vector& x, const nsblowup::Vector& y) {
        return nsblowup::Vector{
            0.08 * x[0] * y[0] + 0.02 * x[1] * y[1],
            -0.03 * x[0] * y[1] + 0.05 * x[1] * y[1]};
    };
    const nsblowup::Vector target{0.04, -0.03};
    auto qtarget = Q(target, target);
    nsblowup::Vector d(2, 0.0);
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) d[i] += B2[i][j] * target[j];
        d[i] += qtarget[i];
    }
    const auto nonlinear = nsblowup::solve_small_quadratic_moment_system(B2, d, Q);
    for (std::size_t i = 0; i < target.size(); ++i)
        if (!close(nonlinear[i], target[i], 5e-10)) return fail("Lemma A.2 fixed-point branch recovery failed");

    // Zero quadratic remainder must reduce exactly to the linear solve, without a smallness restriction.
    const nsblowup::BilinearMap zeroQ = [](const nsblowup::Vector& x, const nsblowup::Vector&) {
        return nsblowup::Vector(x.size(), 0.0);
    };
    const nsblowup::Vector large_d{10.0, -7.0};
    const auto linear_via_fixed_point = nsblowup::solve_small_quadratic_moment_system(B2, large_d, zeroQ);
    const auto linear_direct = nsblowup::solve_linear(B2, large_d);
    for (std::size_t i = 0; i < 2; ++i)
        if (!close(linear_via_fixed_point[i], linear_direct[i])) return fail("Q=0 reduction failed");

    std::cout << "Appendix A moment correction checks passed\n";
    return 0;
}
