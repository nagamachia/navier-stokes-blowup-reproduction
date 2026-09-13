#pragma once

#include "moment_corrections.hpp"
#include "outer_profile_stages.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace nsblowup {

struct AppendixA7MomentChange {
    double Cp{};  // pressure increment: int E^2/(2X) dX
    double S{};   // with U=0: -1/2 int E^2 dX
    double I{};   // angular moment: int sqrt(2X) E dX

    Vector vector() const { return {Cp, S, I}; }
};

inline double appendix_a7_compact_bump(double y, double center, double half_width) {
    if (!(half_width > 0.0)) throw std::invalid_argument("bump half-width must be positive");
    const double u = (y - center) / half_width;
    if (std::abs(u) >= 1.0) return 0.0;
    return std::exp(1.0 - 1.0 / (1.0 - u * u));
}

// Proposition A.7: three additive E-bumps in the second reserved A.9 patch.
// Work in normalized x=X/X_* coordinates, so X_*=e_*=1 in this map.
// The unmodified intermediate power law is E0=f(eta)x^{-1/2-lambda}, U=0.
class AppendixA7HeatCompensation {
public:
    AppendixA7HeatCompensation(double eta, double lambda, std::size_t panels = 4096)
        : eta_(eta), lambda_(lambda), panels_(panels) {
        if (std::abs(eta_) > 1.0) throw std::invalid_argument("|eta| must not exceed 1");
        if (!(lambda_ > 0.0 && lambda_ < 0.25))
            throw std::invalid_argument("lambda must lie in (0,1/4)");
    }

    double base_e(double x) const {
        if (!(x > 0.0)) throw std::invalid_argument("x must be positive");
        return outer_shape_f(eta_) * std::pow(x, -0.5 - lambda_);
    }

    std::array<std::function<double(double)>, 3> bumps() const {
        // The heat reserved patch has logarithmic width 5. Keep all supports
        // strictly inside it and mutually separated, as required by Lemma A.1.
        return {
            [](double x) { return appendix_a7_compact_bump(std::log(x), 0.8, 0.45); },
            [](double x) { return appendix_a7_compact_bump(std::log(x), 2.5, 0.45); },
            [](double x) { return appendix_a7_compact_bump(std::log(x), 4.2, 0.45); }
        };
    }

    AppendixA7MomentChange change(const Vector& c) const {
        if (c.size() != 3) throw std::invalid_argument("three E coefficients required");
        const auto beta = bumps();
        auto delta_e = [&](double x) {
            return c[0] * beta[0](x) + c[1] * beta[1](x) + c[2] * beta[2](x);
        };
        const double a = 1.0;
        const double b = std::exp(5.0);
        AppendixA7MomentChange out;
        out.Cp = composite_simpson([&](double x) {
            const double e0 = base_e(x), de = delta_e(x);
            return (2.0 * e0 * de + de * de) / (2.0 * x);
        }, a, b, panels_);
        out.S = composite_simpson([&](double x) {
            const double e0 = base_e(x), de = delta_e(x);
            return -0.5 * (2.0 * e0 * de + de * de);
        }, a, b, panels_);
        out.I = composite_simpson([&](double x) {
            return std::sqrt(2.0 * x) * delta_e(x);
        }, a, b, panels_);
        return out;
    }

    Matrix jacobian_at_zero() const {
        const auto beta = bumps();
        Matrix B(3, Vector(3, 0.0));
        const double a = 1.0, b = std::exp(5.0);
        for (std::size_t j = 0; j < 3; ++j) {
            B[0][j] = composite_simpson([&](double x) {
                return base_e(x) * beta[j](x) / x;
            }, a, b, panels_);
            B[1][j] = composite_simpson([&](double x) {
                return -base_e(x) * beta[j](x);
            }, a, b, panels_);
            B[2][j] = composite_simpson([&](double x) {
                return std::sqrt(2.0 * x) * beta[j](x);
            }, a, b, panels_);
        }
        return B;
    }

    Vector quadratic_remainder(const Vector& x, const Vector& y) const {
        if (x.size() != 3 || y.size() != 3) throw std::invalid_argument("three coefficients required");
        const auto beta = bumps();
        const double a = 1.0, b = std::exp(5.0);
        Vector q(3, 0.0);
        q[0] = composite_simpson([&](double r) {
            double dx = 0.0, dy = 0.0;
            for (std::size_t j = 0; j < 3; ++j) {
                dx += x[j] * beta[j](r);
                dy += y[j] * beta[j](r);
            }
            return dx * dy / (2.0 * r);
        }, a, b, panels_);
        q[1] = composite_simpson([&](double r) {
            double dx = 0.0, dy = 0.0;
            for (std::size_t j = 0; j < 3; ++j) {
                dx += x[j] * beta[j](r);
                dy += y[j] * beta[j](r);
            }
            return -0.5 * dx * dy;
        }, a, b, panels_);
        // I is exactly linear.
        return q;
    }

    Vector solve_for_change(const Vector& target, std::size_t max_iterations = 100) const {
        if (target.size() != 3) throw std::invalid_argument("three target moments required");
        const Matrix B = jacobian_at_zero();
        Vector scale(3, 0.0), scaled_target = target;
        Matrix scaled_B = B;
        for (std::size_t i = 0; i < 3; ++i) {
            for (double v : B[i]) scale[i] = std::max(scale[i], std::abs(v));
            if (!(scale[i] > 0.0)) throw std::runtime_error("zero row in A.7 Jacobian");
            for (double& v : scaled_B[i]) v /= scale[i];
            scaled_target[i] /= scale[i];
        }
        BilinearMap Q = [&](const Vector& x, const Vector& y) {
            Vector q = quadratic_remainder(x, y);
            for (std::size_t i = 0; i < 3; ++i) q[i] /= scale[i];
            return q;
        };
        return solve_small_quadratic_moment_system(
            scaled_B, scaled_target, Q, max_iterations, 1e-12);
    }

private:
    double eta_;
    double lambda_;
    std::size_t panels_;
};

}  // namespace nsblowup
