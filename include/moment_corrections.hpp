#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nsblowup {

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
using BilinearMap = std::function<Vector(const Vector&, const Vector&)>;

inline double composite_simpson(
    const std::function<double(double)>& f,
    double a,
    double b,
    std::size_t panels = 1024) {
    if (!(b > a)) throw std::invalid_argument("composite_simpson requires b > a");
    if (panels < 2) panels = 2;
    if (panels % 2 != 0) ++panels;
    const double h = (b - a) / static_cast<double>(panels);
    double sum = f(a) + f(b);
    for (std::size_t i = 1; i < panels; ++i) {
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(a + h * static_cast<double>(i));
    }
    return sum * h / 3.0;
}

// Lemma A.1: B_ij = integral x^{alpha_i} beta_j(x) dx.
// The caller supplies separated compact support intervals I_j=[a_j,b_j]
// and the corresponding nonnegative bump functions beta_j.
inline Matrix power_moment_matrix(
    const Vector& alphas,
    const std::vector<std::pair<double, double>>& intervals,
    const std::vector<std::function<double(double)>>& bumps,
    std::size_t panels = 1024) {
    const std::size_t m = alphas.size();
    if (m == 0 || intervals.size() != m || bumps.size() != m) {
        throw std::invalid_argument("moment matrix inputs must have equal nonzero size");
    }
    for (std::size_t j = 0; j < m; ++j) {
        if (!(intervals[j].first > 0.0 && intervals[j].second > intervals[j].first)) {
            throw std::invalid_argument("moment intervals must lie in (0, infinity)");
        }
        if (j > 0 && !(intervals[j - 1].second < intervals[j].first)) {
            throw std::invalid_argument("moment intervals must be strictly separated and ordered");
        }
    }

    Matrix B(m, Vector(m, 0.0));
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            const auto [a, b] = intervals[j];
            B[i][j] = composite_simpson(
                [&](double x) { return std::pow(x, alphas[i]) * bumps[j](x); },
                a, b, panels);
        }
    }
    return B;
}

inline Vector solve_linear(Matrix A, Vector b, double pivot_tolerance = 1e-13) {
    const std::size_t n = A.size();
    if (n == 0 || b.size() != n) throw std::invalid_argument("invalid linear system size");
    for (const auto& row : A) if (row.size() != n) throw std::invalid_argument("matrix must be square");

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        }
        if (std::abs(A[pivot][k]) <= pivot_tolerance) {
            throw std::runtime_error("moment Jacobian is numerically singular");
        }
        if (pivot != k) {
            std::swap(A[pivot], A[k]);
            std::swap(b[pivot], b[k]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = A[i][k] / A[k][k];
            A[i][k] = 0.0;
            for (std::size_t j = k + 1; j < n; ++j) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    Vector x(n, 0.0);
    for (std::size_t ii = n; ii-- > 0;) {
        double rhs = b[ii];
        for (std::size_t j = ii + 1; j < n; ++j) rhs -= A[ii][j] * x[j];
        x[ii] = rhs / A[ii][ii];
    }
    return x;
}

inline double norm_inf(const Vector& x) {
    double result = 0.0;
    for (double v : x) result = std::max(result, std::abs(v));
    return result;
}

// Lemma A.2 fixed-point branch:
// c <- B^{-1}(d - Q(c,c)).
// This routine deliberately starts from zero, matching the branch selected in the paper.
inline Vector solve_small_quadratic_moment_system(
    const Matrix& B,
    const Vector& discrepancy,
    const BilinearMap& Q,
    std::size_t max_iterations = 100,
    double tolerance = 1e-12) {
    if (B.empty() || discrepancy.size() != B.size()) {
        throw std::invalid_argument("quadratic moment system has inconsistent dimensions");
    }
    Vector c(B.size(), 0.0);
    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
        const Vector q = Q(c, c);
        if (q.size() != c.size()) throw std::invalid_argument("bilinear map returned wrong dimension");
        Vector rhs = discrepancy;
        for (std::size_t i = 0; i < rhs.size(); ++i) rhs[i] -= q[i];
        Vector next = solve_linear(B, rhs);
        Vector delta(next.size());
        for (std::size_t i = 0; i < next.size(); ++i) delta[i] = next[i] - c[i];
        c = std::move(next);
        if (norm_inf(delta) <= tolerance * std::max(1.0, norm_inf(c))) return c;
    }
    throw std::runtime_error("quadratic moment correction fixed point did not converge");
}

inline double determinant(Matrix A, double pivot_tolerance = 1e-15) {
    const std::size_t n = A.size();
    if (n == 0) return 1.0;
    for (const auto& row : A) if (row.size() != n) throw std::invalid_argument("matrix must be square");
    double det = 1.0;
    int sign = 1;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i) if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        if (std::abs(A[pivot][k]) <= pivot_tolerance) return 0.0;
        if (pivot != k) { std::swap(A[pivot], A[k]); sign = -sign; }
        const double p = A[k][k];
        det *= p;
        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = A[i][k] / p;
            for (std::size_t j = k + 1; j < n; ++j) A[i][j] -= factor * A[k][j];
        }
    }
    return static_cast<double>(sign) * det;
}

}  // namespace nsblowup
