#pragma once

#include "moment_corrections.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace nsblowup {

using ScalarProfile1d = std::function<double(double)>;

struct FiveMomentChange {
    double M{};
    double I{};
    double J{};
    double S{};
    double Cp{};

    Vector vector() const { return {M, I, J, S, Cp}; }
};

// Numerical realization of the exact degree-two five-moment map used in
// Corollary A.3. The correction patch is strictly away from X=0.
class FiveBumpMomentPatch {
public:
    FiveBumpMomentPatch(
        double a,
        double b,
        ScalarProfile1d base_u,
        ScalarProfile1d base_e,
        std::array<ScalarProfile1d, 2> u_bumps,
        std::array<ScalarProfile1d, 3> e_bumps,
        std::size_t panels = 2048)
        : a_(a), b_(b), base_u_(std::move(base_u)), base_e_(std::move(base_e)),
          u_bumps_(std::move(u_bumps)), e_bumps_(std::move(e_bumps)), panels_(panels) {
        if (!(a_ > 0.0 && b_ > a_)) throw std::invalid_argument("correction patch must lie in X>0");
    }

    FiveMomentChange change(const Vector& c) const {
        if (c.size() != 5) throw std::invalid_argument("five correction coefficients required");
        FiveMomentChange out;
        out.M = composite_simpson([&](double x) {
            return corrected_u(x, c) - base_u_(x);
        }, a_, b_, panels_);
        out.I = composite_simpson([&](double x) {
            return std::sqrt(2.0 * x) * (corrected_e(x, c) - base_e_(x));
        }, a_, b_, panels_);
        out.J = composite_simpson([&](double x) {
            const double u = corrected_u(x, c), e = corrected_e(x, c);
            return std::sqrt(2.0 * x) * (u * e - base_u_(x) * base_e_(x));
        }, a_, b_, panels_);
        out.S = composite_simpson([&](double x) {
            const double u = corrected_u(x, c), e = corrected_e(x, c);
            const double u0 = base_u_(x), e0 = base_e_(x);
            return u * u - u0 * u0 - 0.5 * (e * e - e0 * e0);
        }, a_, b_, panels_);
        out.Cp = composite_simpson([&](double x) {
            const double e = corrected_e(x, c), e0 = base_e_(x);
            return (e * e - e0 * e0) / (2.0 * x);
        }, a_, b_, panels_);
        return out;
    }

    Matrix jacobian_at_zero() const {
        Matrix B(5, Vector(5, 0.0));
        for (std::size_t col = 0; col < 5; ++col) {
            FiveMomentChange d;
            if (col < 2) {
                const auto& beta = u_bumps_[col];
                d.M = composite_simpson([&](double x) { return beta(x); }, a_, b_, panels_);
                d.J = composite_simpson([&](double x) {
                    return std::sqrt(2.0 * x) * base_e_(x) * beta(x);
                }, a_, b_, panels_);
                d.S = composite_simpson([&](double x) {
                    return 2.0 * base_u_(x) * beta(x);
                }, a_, b_, panels_);
            } else {
                const auto& beta = e_bumps_[col - 2];
                d.I = composite_simpson([&](double x) {
                    return std::sqrt(2.0 * x) * beta(x);
                }, a_, b_, panels_);
                d.J = composite_simpson([&](double x) {
                    return std::sqrt(2.0 * x) * base_u_(x) * beta(x);
                }, a_, b_, panels_);
                d.S = composite_simpson([&](double x) {
                    return -base_e_(x) * beta(x);
                }, a_, b_, panels_);
                d.Cp = composite_simpson([&](double x) {
                    return base_e_(x) * beta(x) / x;
                }, a_, b_, panels_);
            }
            const Vector v = d.vector();
            for (std::size_t row = 0; row < 5; ++row) B[row][col] = v[row];
        }
        return B;
    }

    Vector quadratic_remainder(const Vector& x, const Vector& y) const {
        if (x.size() != 5 || y.size() != 5) throw std::invalid_argument("five coefficients required");
        Vector xy(5), xpy(5);
        for (std::size_t i = 0; i < 5; ++i) xpy[i] = x[i] + y[i];
        const Matrix B = jacobian_at_zero();
        auto remainder = [&](const Vector& c) {
            Vector r = change(c).vector();
            for (std::size_t i = 0; i < 5; ++i)
                for (std::size_t j = 0; j < 5; ++j) r[i] -= B[i][j] * c[j];
            return r;
        };
        const Vector rx = remainder(x), ry = remainder(y), rxy = remainder(xpy);
        for (std::size_t i = 0; i < 5; ++i) xy[i] = 0.5 * (rxy[i] - rx[i] - ry[i]);
        return xy;
    }

    Vector solve_for_change(const Vector& discrepancy, std::size_t max_iterations = 100) const {
        if (discrepancy.size() != 5) throw std::invalid_argument("five moment discrepancies required");
        const Matrix B = jacobian_at_zero();

        // The five rows have different physical scales. Row-normalizing avoids
        // mistaking a small but perfectly invertible Jacobian for a singular one.
        Vector row_scale(5, 0.0);
        Matrix scaled_B = B;
        Vector scaled_d = discrepancy;
        for (std::size_t i = 0; i < 5; ++i) {
            for (double value : B[i]) row_scale[i] = std::max(row_scale[i], std::abs(value));
            if (!(row_scale[i] > 0.0)) throw std::runtime_error("zero row in five-moment Jacobian");
            for (double& value : scaled_B[i]) value /= row_scale[i];
            scaled_d[i] /= row_scale[i];
        }

        BilinearMap scaled_Q = [&](const Vector& x, const Vector& y) {
            Vector q = quadratic_remainder(x, y);
            for (std::size_t i = 0; i < q.size(); ++i) q[i] /= row_scale[i];
            return q;
        };
        return solve_small_quadratic_moment_system(scaled_B, scaled_d, scaled_Q, max_iterations, 1e-11);
    }

private:
    double corrected_u(double x, const Vector& c) const {
        return base_u_(x) + c[0] * u_bumps_[0](x) + c[1] * u_bumps_[1](x);
    }
    double corrected_e(double x, const Vector& c) const {
        return base_e_(x) + c[2] * e_bumps_[0](x) + c[3] * e_bumps_[1](x) + c[4] * e_bumps_[2](x);
    }

    double a_, b_;
    ScalarProfile1d base_u_, base_e_;
    std::array<ScalarProfile1d, 2> u_bumps_;
    std::array<ScalarProfile1d, 3> e_bumps_;
    std::size_t panels_;
};

}  // namespace nsblowup
