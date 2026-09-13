#pragma once

#include "appendix_b_axis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

struct AppendixBEndpointAudit {
    double min_phi0{};
    double min_azimuthal_p1_on_chi_branch{};
    double min_abs_normalized_ns_on_axial_branch{};
    double max_log_required_C_normalized{};
    bool azimuthal_branch_certified{};
    bool axial_branch_separated{};
};

// f0(z)+z f0'(z)=sum (-z/2)^a/(a!)^2.
inline double appendix_b_f0_plus_z_df0(double z, int terms = 80) {
    if (!(z >= 0.0)) throw std::invalid_argument("Appendix B requires z>=0");
    double sum = 1.0;
    double term = 1.0;
    for (int a = 1; a < terms; ++a) {
        term *= (-0.5 * z) / (static_cast<double>(a) * static_cast<double>(a));
        sum += term;
        if (std::abs(term) < 1e-16 * std::max(1.0, std::abs(sum))) break;
    }
    return sum;
}

inline double appendix_b_unperturbed_p1(double Y, double chi) {
    const double z = Y * chi;
    const double f = appendix_b_f0(z);
    const double fplus = appendix_b_f0_plus_z_df0(z);
    return 2.0 - 2.0 * fplus / f;
}

// Proposition B.3 endpoint alternative in the Lambda->infinity comparison model.
// On chi>.99, the azimuthal term alone gives p1>2.36 at Y=4.
// On chi<=.99, B.2 gives |Z*|>delta*. We record the finite amplitude scale
// required for p2^2/p1>2.3 in normalized pressure units. This is an asymptotic
// audit of the mechanism, not a replacement for the full nonlinear B.2 solve.
inline AppendixBEndpointAudit appendix_b_endpoint_asymptotic_audit(
    const OuterProfileParameters& p,
    const AppendixBSeparationAudit& sep,
    double Lambda = 25.0, double j0 = 0.05,
    int eta_samples = 241) {
    if (!(Lambda >= 1.0)) throw std::invalid_argument("Appendix B requires Lambda>=1");
    const double h = std::exp(p.log_h);
    double min_phi0 = 1.0;
    double min_p1 = std::numeric_limits<double>::infinity();
    double min_abs_ns = std::numeric_limits<double>::infinity();
    double max_logC = -std::numeric_limits<double>::infinity();
    bool saw_chi = false, saw_axial = false;

    for (int i = 0; i < eta_samples; ++i) {
        const double eta = -1.0 + 2.0 * static_cast<double>(i) /
                                      static_cast<double>(eta_samples - 1);
        const double H = appendix_b_Hstar(eta, h, j0);
        const double L = appendix_b_L(eta, h);
        const double chi = H * H / (H * H + sep.sigma_star * sep.sigma_star);
        const double phi0 = appendix_b_f0(4.0 * chi);
        min_phi0 = std::min(min_phi0, phi0);

        if (chi > 0.99) {
            saw_chi = true;
            min_p1 = std::min(min_p1, appendix_b_unperturbed_p1(4.0, chi));
        } else {
            const auto a = appendix_b_axis_data(p, eta, j0, 20.0, 0.05, 0.04);
            if (std::abs(a.normalized_Zstar) <= sep.normalized_delta_star)
                throw std::runtime_error("B.2 separation inconsistent with endpoint branch");
            saw_axial = true;
            const double ns_norm = std::abs(a.normalized_Zstar) / L;
            min_abs_ns = std::min(min_abs_ns, ns_norm);
            const double p1 = std::max(appendix_b_unperturbed_p1(4.0, chi), 1e-12);
            const double log_phi = appendix_b_log_phi_star(
                p, eta, Lambda, sep.sigma_star, j0, 400);
            // |p2| = C sqrt(2/Lambda) |ns|/(phi*Phi).  Z and ns are
            // normalized by P_*^2 here, so C is reported in the same normalized audit.
            const double log_required = 0.5 * std::log(2.3 * p1)
                - 0.5 * std::log(2.0 / Lambda)
                - std::log(ns_norm) + log_phi + std::log(phi0);
            max_logC = std::max(max_logC, log_required);
        }
    }

    AppendixBEndpointAudit out;
    out.min_phi0 = min_phi0;
    out.min_azimuthal_p1_on_chi_branch = saw_chi ? min_p1 : 0.0;
    out.min_abs_normalized_ns_on_axial_branch = saw_axial ? min_abs_ns : 0.0;
    out.max_log_required_C_normalized = saw_axial ? max_logC : 0.0;
    out.azimuthal_branch_certified = saw_chi && min_p1 > 2.36;
    out.axial_branch_separated = saw_axial && min_abs_ns > sep.normalized_delta_star;
    return out;
}

} // namespace nsblowup
