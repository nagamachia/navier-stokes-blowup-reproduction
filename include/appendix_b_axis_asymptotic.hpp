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
    double min_abs_normalized_Z_on_axial_branch{};
    double min_abs_normalized_ns_on_axial_branch{};
    double max_log_required_C_normalized{};
    bool azimuthal_branch_certified{};
    bool axial_branch_separated{};
};

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

inline AppendixBEndpointAudit appendix_b_endpoint_asymptotic_audit(
    const OuterProfileParameters& p,
    const AppendixBSeparationAudit& sep,
    double Lambda = 25.0, double j0 = 0.05,
    int eta_samples = 241) {
    if (!(Lambda >= 1.0)) throw std::invalid_argument("Appendix B requires Lambda>=1");
    const double h = std::exp(p.log_h);
    double min_phi0 = 1.0;
    double min_p1 = std::numeric_limits<double>::infinity();
    double min_abs_Z = std::numeric_limits<double>::infinity();
    double min_abs_ns = std::numeric_limits<double>::infinity();
    double max_logC = -std::numeric_limits<double>::infinity();
    bool saw_chi = false, saw_axial = false;

    auto audit_eta = [&](double eta) {
        const double H = appendix_b_Hstar(eta, h, j0);
        const double L = appendix_b_L(eta, h);
        const double chi = H * H / (H * H + sep.sigma_star * sep.sigma_star);
        const double phi0 = appendix_b_f0(4.0 * chi);
        min_phi0 = std::min(min_phi0, phi0);
        if (chi > 0.99) {
            saw_chi = true;
            min_p1 = std::min(min_p1, appendix_b_unperturbed_p1(4.0, chi));
            return;
        }

        const auto a = appendix_b_axis_data(p, eta, j0, 20.0, 0.05, 0.04);
        const double abs_Z = std::abs(a.normalized_Zstar);
        if (!(abs_Z > sep.normalized_delta_star))
            throw std::runtime_error("B.2 separation inconsistent with endpoint branch");
        saw_axial = true;
        min_abs_Z = std::min(min_abs_Z, abs_Z);
        const double ns_norm = abs_Z / L;
        min_abs_ns = std::min(min_abs_ns, ns_norm);
        const double p1 = std::max(appendix_b_unperturbed_p1(4.0, chi), 1e-12);
        const double log_phi = appendix_b_log_phi_star(
            p, eta, Lambda, sep.sigma_star, j0, 400);
        const double log_required = 0.5 * std::log(2.3 * p1)
            - 0.5 * std::log(2.0 / Lambda)
            - std::log(ns_norm) + log_phi + std::log(phi0);
        max_logC = std::max(max_logC, log_required);
    };

    for (int i = 0; i < eta_samples; ++i) {
        const double eta = -1.0 + 2.0 * static_cast<double>(i) /
                                      static_cast<double>(eta_samples - 1);
        audit_eta(eta);
    }
    // The chi<=.99 set can be much narrower than a uniform eta grid because
    // sigma_* is deliberately tiny. Include the unique zero H*(eta0)=0 exactly;
    // here chi=0, so this point belongs to the axial alternative by construction.
    audit_eta(sep.eta0);

    AppendixBEndpointAudit out;
    out.min_phi0 = min_phi0;
    out.min_azimuthal_p1_on_chi_branch = saw_chi ? min_p1 : 0.0;
    out.min_abs_normalized_Z_on_axial_branch = saw_axial ? min_abs_Z : 0.0;
    out.min_abs_normalized_ns_on_axial_branch = saw_axial ? min_abs_ns : 0.0;
    out.max_log_required_C_normalized = saw_axial ? max_logC : 0.0;
    out.azimuthal_branch_certified = saw_chi && min_p1 > 2.36;
    out.axial_branch_separated = saw_axial && min_abs_Z > sep.normalized_delta_star;
    return out;
}

} // namespace nsblowup
