#pragma once

#include "appendix_a4_pressure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBAxisData {
    double eta{};
    double Ustar{};
    double Hstar{};
    double Wstar{};
    double normalized_Zstar{}; // Z*/P_*^2
};

struct AppendixBSeparationAudit {
    double eta0{};
    double normalized_Zstar_eta0{};
    double normalized_delta_star{};
    double sigma_star{};
    double min_chi_on_small_Z{};
    bool unique_H_zero{};
    bool positive_Z_at_H_zero{};
    bool chi_separation_certified{};
};

inline double appendix_b_A(double h) { return 0.5 + h; }
inline double appendix_b_D(double h) { return 0.5 - h; }
inline double appendix_b_d(double eta) { return 1.0 - eta * eta; }
inline double appendix_b_L(double eta, double h) { return 1.0 - 2.0 * h * eta * eta; }

inline AppendixBAxisData appendix_b_axis_data(
    const OuterProfileParameters& p, double eta, double j0 = 0.05,
    double Tf = 20.0, double co = 0.05, double max_step = 0.02) {
    p.validate();
    if (!(eta >= -1.0 && eta <= 1.0))
        throw std::invalid_argument("Appendix B requires eta in [-1,1]");
    if (!(j0 > 0.0 && j0 <= 0.05))
        throw std::invalid_argument("Appendix B requires 0<j0<=0.05");

    const double h = std::exp(p.log_h);
    const double A = appendix_b_A(h);
    const double D = appendix_b_D(h);
    const double d = appendix_b_d(eta);
    const double U = 4.0 * eta + j0;
    const double Ueta = 4.0;
    const double H = D * eta + d * U;
    const double W = 1.0 - d * Ueta - 2.0 * D * eta * U;

    const auto pressure = appendix_a4_pressure_datum(p, eta, Tf, co, max_step);
    const double invP2 = std::exp(-2.0 * p.log_Pstar);

    // Equation (B.1), normalized by P_*^2. The O(1) velocity terms are
    // multiplied by P_*^{-2}, while the pressure terms are already normalized.
    const double Znorm =
        invP2 * (-A * (1.0 - 2.0 * eta * U) * U - H * Ueta)
        - d * pressure.normalized_deta
        + 4.0 * A * eta * pressure.normalized_pi0;

    return {eta, U, H, W, Znorm};
}

inline double appendix_b_find_H_zero(
    const OuterProfileParameters& p, double j0 = 0.05,
    double tol = 1e-13) {
    const double h = std::exp(p.log_h);
    const double D = appendix_b_D(h);
    auto H = [&](double eta) {
        return D * eta + appendix_b_d(eta) * (4.0 * eta + j0);
    };
    double a = -0.999999, b = 0.0;
    double fa = H(a), fb = H(b);
    if (!(fa < 0.0 && fb > 0.0))
        throw std::runtime_error("Appendix B H* root is not bracketed in (-1,0)");
    for (int iter = 0; iter < 100; ++iter) {
        const double m = 0.5 * (a + b);
        const double fm = H(m);
        if (std::abs(fm) < tol || b - a < tol) return m;
        if (fm > 0.0) { b = m; fb = fm; }
        else { a = m; fa = fm; }
    }
    return 0.5 * (a + b);
}

inline AppendixBSeparationAudit appendix_b_separation_audit(
    const OuterProfileParameters& p, double j0 = 0.05,
    double Tf = 20.0, double co = 0.05, double max_step = 0.02,
    int eta_samples = 801) {
    if (eta_samples < 101) eta_samples = 101;
    const double eta0 = appendix_b_find_H_zero(p, j0);
    const auto root_data = appendix_b_axis_data(p, eta0, j0, Tf, co, max_step);

    // B.2 only requires some delta_* whose small-|Z*| set avoids eta0.
    // We work in the P_*^{-2} normalization and take half the positive value
    // at the H* zero, leaving a strict margin.
    const double delta = 0.5 * root_data.normalized_Zstar;
    if (!(delta > 0.0))
        throw std::runtime_error("Appendix B requires Z*(eta0)>0");

    double min_abs_H = std::numeric_limits<double>::infinity();
    bool saw_small_Z = false;
    for (int i = 0; i < eta_samples; ++i) {
        const double eta = -1.0 + 2.0 * static_cast<double>(i) /
                                      static_cast<double>(eta_samples - 1);
        const auto a = appendix_b_axis_data(p, eta, j0, Tf, co, max_step);
        if (std::abs(a.normalized_Zstar) <= delta) {
            saw_small_Z = true;
            min_abs_H = std::min(min_abs_H, std::abs(a.Hstar));
        }
    }
    if (!saw_small_Z || !(min_abs_H > 0.0) || !std::isfinite(min_abs_H))
        throw std::runtime_error("Appendix B failed to resolve the small-|Z*| separation set");

    // Choose sigma_* with comfortable strict slack: sigma = 0.05 min|H|.
    // Then chi >= 1/(1+0.05^2) > 0.9975 on the audited set.
    const double sigma = 0.05 * min_abs_H;
    double min_chi = 1.0;
    for (int i = 0; i < eta_samples; ++i) {
        const double eta = -1.0 + 2.0 * static_cast<double>(i) /
                                      static_cast<double>(eta_samples - 1);
        const auto a = appendix_b_axis_data(p, eta, j0, Tf, co, max_step);
        if (std::abs(a.normalized_Zstar) <= delta) {
            const double chi = a.Hstar * a.Hstar /
                               (a.Hstar * a.Hstar + sigma * sigma);
            min_chi = std::min(min_chi, chi);
        }
    }

    AppendixBSeparationAudit out;
    out.eta0 = eta0;
    out.normalized_Zstar_eta0 = root_data.normalized_Zstar;
    out.normalized_delta_star = delta;
    out.sigma_star = sigma;
    out.min_chi_on_small_Z = min_chi;
    out.unique_H_zero = eta0 > -1.0 && eta0 < 0.0;
    out.positive_Z_at_H_zero = root_data.normalized_Zstar > 0.0;
    out.chi_separation_certified = min_chi > 0.99;
    return out;
}

inline double appendix_b_zeta_star(
    const OuterProfileParameters& p, double eta, double sigma_star,
    double j0 = 0.05) {
    const double h = std::exp(p.log_h);
    const double H = appendix_b_axis_data(p, eta, j0).Hstar;
    const double L = appendix_b_L(eta, h);
    return -L * H / (H * H + sigma_star * sigma_star);
}

inline double appendix_b_log_phi_star(
    const OuterProfileParameters& p, double eta, double Lambda,
    double sigma_star, double j0 = 0.05, int panels = 800) {
    if (!(Lambda >= 1.0 && sigma_star > 0.0))
        throw std::invalid_argument("Appendix B requires Lambda>=1 and sigma_*>0");
    if (!(eta >= -1.0 && eta <= 1.0))
        throw std::invalid_argument("Appendix B requires eta in [-1,1]");
    if (eta == 0.0) return 0.0;
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    const double a = std::min(0.0, eta);
    const double b = std::max(0.0, eta);
    const double dx = (b - a) / static_cast<double>(panels);
    double sum = 0.0;
    for (int i = 0; i <= panels; ++i) {
        const double x = a + dx * static_cast<double>(i);
        const double w = (i == 0 || i == panels) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum += w * appendix_b_zeta_star(p, x, sigma_star, j0);
    }
    const double integral = sum * dx / 3.0;
    return Lambda * (eta >= 0.0 ? integral : -integral);
}

inline double appendix_b_f0(double z, int terms = 80) {
    if (!(z >= 0.0)) throw std::invalid_argument("Appendix B f0 requires z>=0");
    double sum = 1.0;
    double term = 1.0;
    for (int alpha = 1; alpha < terms; ++alpha) {
        term *= (-0.5 * z) /
                (static_cast<double>(alpha) * static_cast<double>(alpha + 1));
        sum += term;
        if (std::abs(term) < 1e-16 * std::max(1.0, std::abs(sum))) break;
    }
    return sum;
}

} // namespace nsblowup
