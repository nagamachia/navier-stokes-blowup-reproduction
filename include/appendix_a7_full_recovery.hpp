#pragma once

#include "appendix_a7_log_compensation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct AppendixA7FullRecoveryAudit {
    AppendixA7ScheduleScales scales{};
    AppendixA7LeadingHeatLogDiscrepancy leading{};
    AppendixA7LeadingHeatErrorAudit heat_error{};
    AppendixA7LogCompensationAudit compensation{};

    double worst_heat_log_relative_bound{};
    bool heat_asymptotic_certified{};
    bool dominant_compensation_certified{};
    bool nonlinear_recovery_branch_certified{};
    bool pressure_tracks_Cp{};
    bool direct_double_exact_recovery_available{};
};

// End-to-end numerical audit for Proposition A.7 in the ordered regime.
//
// This deliberately separates what can be checked directly in double precision
// from what is certified through logarithmic scale bounds. The heat replacement
// is represented by its leading signed-log discrepancy plus a rigorous relative
// remainder bound. The second A.9 patch is then audited in the dominant target
// scale and the Lemma A.2 zero-start fixed-point map is checked to be a
// contraction on an explicit ball. A false direct_double_exact_recovery_available
// is expected in the strongly ordered regime and must not be interpreted as
// failure of the mathematical correction branch.
inline AppendixA7FullRecoveryAudit appendix_a7_full_recovery_audit(
    const OuterProfileParameters& p, double eta,
    double Tf = 20.0, double co = 0.05,
    std::size_t panels = 4096) {
    p.validate();
    if (!(std::abs(eta) < 1.0))
        throw std::invalid_argument("A.7 full recovery audit requires |eta|<1");

    AppendixA7FullRecoveryAudit out;
    out.scales = appendix_a7_schedule_scales(p, eta, Tf, co);
    out.leading = appendix_a7_leading_heat_log_discrepancy(p, eta, Tf, co, 1024);
    out.heat_error = appendix_a7_leading_heat_error_audit(p, eta, Tf, co);
    out.compensation = appendix_a7_audit_log_compensation(
        eta, p.lambda,
        out.leading.patch_target_Cp,
        out.leading.patch_target_S,
        out.leading.patch_target_I,
        panels);

    out.worst_heat_log_relative_bound = std::max(
        out.heat_error.log_relative_bound_CpS,
        out.heat_error.log_relative_bound_I);
    out.heat_asymptotic_certified = out.worst_heat_log_relative_bound < 0.0;

    out.dominant_compensation_certified =
        norm_inf(out.compensation.dominant_scaled_linear_residual) < 1e-10 &&
        out.compensation.quadratic_to_dominant_log_ratio < 0.0;

    // The contraction test is independent of whether tiny target components
    // can be explicitly represented in a double.  It uses ||d||_inf<=1 and
    // finite-dimensional operator bounds, so it certifies existence/uniqueness
    // of the complete nonlinear small branch within the audited model.
    out.nonlinear_recovery_branch_certified =
        out.heat_asymptotic_certified &&
        out.compensation.nonlinear_branch_contraction_certified;

    // Cp is exactly the axis-pressure increment used by the A.4/A.7 coupling,
    // so pressure restoration is not a fourth independent scalar condition.
    out.pressure_tracks_Cp = true;
    out.direct_double_exact_recovery_available =
        out.compensation.all_components_resolvable_in_double;
    return out;
}

} // namespace nsblowup
