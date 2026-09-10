#pragma once

#include "moment_corrections.hpp"
#include "outer_profile_stages.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace nsblowup {

inline double appendix_a3_simpson(const std::function<double(double)>& f,
                                  double a,
                                  double b,
                                  int panels = 1024) {
    if (a == b) return 0.0;
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    const double h = (b - a) / static_cast<double>(panels);
    double sum = f(a) + f(b);
    for (int i = 1; i < panels; ++i) {
        const double x = a + h * static_cast<double>(i);
        sum += (i % 2 ? 4.0 : 2.0) * f(x);
    }
    return sum * h / 3.0;
}

// Section A.2/A.3 axial pulse. Because sigma(v/.02)=1 for v>=.02 and
// sigma(s)+sigma(1-s)=1, integral_0^.02 sigma(v/.02)dv=.01 exactly.
inline double appendix_a3_phi_b(double xi) {
    if (xi <= 0.0) return 0.0;
    if (xi >= 0.02) return xi - 0.01;
    return 0.02 * appendix_a3_simpson(
        [](double s) { return appendix_a_smooth_step(s); },
        0.0, xi / 0.02, 200);
}

inline double appendix_a3_R0(double xi) {
    if (xi <= 0.0 || xi >= 11.0) return 0.0;
    return appendix_a3_phi_b(xi) * (1.0 - appendix_a_smooth_step(xi - 10.0));
}

// K_b in (A.19).
inline double appendix_a3_Kb(int panels = 6000) {
    return appendix_a3_simpson(
        [](double xi) {
            const double r = appendix_a3_R0(xi);
            return std::exp(-2.0 * xi) * r * r;
        },
        0.0, 13.0, panels);
}

// Identical translates of a width-w rescaled sigma' bump; integral beta dy = 1.
inline double appendix_a3_unit_bump(double y, double center, double width = 0.3) {
    if (!(width > 0.0)) throw std::invalid_argument("bump width must be positive");
    const double left = center - 0.5 * width;
    const double s = (y - left) / width;
    if (s <= 0.0 || s >= 1.0) return 0.0;
    return appendix_a_smooth_step_derivative(s) / width;
}

struct PulseMomentClosure {
    double c1{};
    double c2{};
    double residual_M{};
    double residual_J{};
};

// Numerical form of (A.15). Inputs pre_M and pre_J are the normalized moment
// discrepancies at the first correction-bump center, after including the main
// pulse with the selected amplitude. The two late bumps are then solved exactly
// in the linear M,J rows. The weights are exp(s_i(y-y1)), s1=1/2-lambda,
// s2=1/2-2lambda, as in the proof of A.15.
inline PulseMomentClosure appendix_a3_close_pulse_MJ(double lambda,
                                                     double pre_M,
                                                     double pre_J,
                                                     double amplitude) {
    if (!(lambda > 0.0 && lambda < 0.2))
        throw std::invalid_argument("lambda must be small and positive");
    const double y1 = 13.0 / lambda - 3.0;
    const double y2 = 13.0 / lambda - 1.0;
    const double s1 = 0.5 - lambda;
    const double s2 = 0.5 - 2.0 * lambda;

    auto main_moment = [&](double slope) {
        // Change xi=lambda*y. The normalized exponential is evaluated relative
        // to y1, avoiding the theorem-scale overflow discussed in Appendix A.
        return appendix_a3_simpson(
            [&](double xi) {
                const double y = xi / lambda;
                return std::exp(slope * (y - y1)) * appendix_a3_R0(xi) / lambda;
            }, 0.0, 11.0, 5000);
    };

    const double rhs_M = -(pre_M + amplitude * main_moment(s1));
    const double rhs_J = -(pre_J + amplitude * main_moment(s2));

    Matrix B(2, Vector(2, 0.0));
    for (int j = 0; j < 2; ++j) {
        const double center = (j == 0 ? y1 : y2);
        const double a = center - 0.15;
        const double b = center + 0.15;
        B[0][j] = appendix_a3_simpson(
            [&](double y) {
                return std::exp(s1 * (y - y1)) * appendix_a3_unit_bump(y, center);
            }, a, b, 800);
        B[1][j] = appendix_a3_simpson(
            [&](double y) {
                return std::exp(s2 * (y - y1)) * appendix_a3_unit_bump(y, center);
            }, a, b, 800);
    }

    const Vector c = solve_linear(B, {rhs_M, rhs_J}, 1e-14);
    PulseMomentClosure out;
    out.c1 = c[0];
    out.c2 = c[1];
    out.residual_M = pre_M + amplitude * main_moment(s1) + B[0][0] * c[0] + B[0][1] * c[1];
    out.residual_J = pre_J + amplitude * main_moment(s2) + B[1][0] * c[0] + B[1][1] * c[1];
    return out;
}

struct AngularMomentClosure {
    double c1{};
    double c2{};
    double residual_I{};
    double residual_pressure{};
};

// Specialized numerical form of the two relative-E corrections in (A.11).
// dI and dP are normalized desired changes at the first bump center. The I row
// is linear. The pressure row is exactly quadratic in the relative correction.
inline AngularMomentClosure appendix_a3_close_angular_I_pressure(double lambda,
                                                                 double dI,
                                                                 double dP) {
    if (!(lambda > 0.0 && lambda < 0.2))
        throw std::invalid_argument("lambda must be small and positive");
    const double length = 30.0 * std::log(1.0 / lambda);
    const double y1 = length - 3.0;
    const double y2 = length - 1.0;
    const double sI = 1.0 - lambda;
    const double sP = -1.0 - 2.0 * lambda;

    auto beta = [&](int j, double y) {
        return appendix_a3_unit_bump(y, j == 0 ? y1 : y2);
    };

    Matrix B(2, Vector(2, 0.0));
    for (int j = 0; j < 2; ++j) {
        const double center = (j == 0 ? y1 : y2);
        const double a = center - 0.15;
        const double b = center + 0.15;
        B[0][j] = appendix_a3_simpson(
            [&](double y) { return std::exp(sI * (y - y1)) * beta(j, y); },
            a, b, 800);
        B[1][j] = appendix_a3_simpson(
            [&](double y) { return 2.0 * std::exp(sP * (y - y1)) * beta(j, y); },
            a, b, 800);
    }

    BilinearMap Q = [&](const Vector& x, const Vector& y) {
        Vector out(2, 0.0);
        const double a = y1 - 0.15;
        const double b = y2 + 0.15;
        out[1] = appendix_a3_simpson(
            [&](double z) {
                const double bx = x[0] * beta(0, z) + x[1] * beta(1, z);
                const double by = y[0] * beta(0, z) + y[1] * beta(1, z);
                return std::exp(sP * (z - y1)) * bx * by;
            }, a, b, 2400);
        return out;
    };

    const Vector c = solve_small_quadratic_moment_system(B, {dI, dP}, Q, 100, 1e-12);
    const Vector q = Q(c, c);
    AngularMomentClosure out;
    out.c1 = c[0];
    out.c2 = c[1];
    out.residual_I = B[0][0] * c[0] + B[0][1] * c[1] - dI;
    out.residual_pressure = B[1][0] * c[0] + B[1][1] * c[1] + q[1] - dP;
    return out;
}

using A3AmplitudeError = std::function<double(double, double)>;

// Solve (A.19) on the paper's bracket [0.9,1.2]. The callback supplies the
// small remainder E(Amp,eta) after dividing by the normalization in (A.19).
inline double appendix_a3_solve_amplitude(double eta,
                                          const A3AmplitudeError& error = {},
                                          int iterations = 80) {
    const double Kb = appendix_a3_Kb();
    const double constant = (1.0 - std::exp(-26.0)) / 4.0;
    auto f = [&](double amp) {
        return amp * amp * Kb - constant + (error ? error(amp, eta) : 0.0);
    };
    double lo = 0.9, hi = 1.2;
    double flo = f(lo), fhi = f(hi);
    if (!(flo < 0.0 && fhi > 0.0))
        throw std::runtime_error("A.19 amplitude root is not bracketed on [0.9,1.2]");
    for (int i = 0; i < iterations; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double fm = f(mid);
        if (fm > 0.0) hi = mid;
        else lo = mid;
    }
    return 0.5 * (lo + hi);
}

// Exact terminal representation (A.16). Q_s(0)=Q_p and Q_s(3)=0.
inline double appendix_a3_terminal_Qs(double y, double h, double rho_o, int panels = 800) {
    if (y <= 0.0) y = 0.0;
    if (y >= 3.0) return 0.0;
    auto exponent = [&](double v) {
        return appendix_a3_simpson(
            [&](double s) { return 1.0 + outer_terminal_l(s, h, rho_o); },
            y, v, 300);
    };
    return appendix_a3_simpson(
        [&](double v) {
            return std::exp(exponent(v)) * outer_terminal_fprime(v, rho_o) /
                   outer_terminal_f(v, rho_o);
        }, y, 3.0, panels);
}

}  // namespace nsblowup
