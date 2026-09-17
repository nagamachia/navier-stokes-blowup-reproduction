#pragma once

#include "appendix_b_Z_interval.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nsblowup {

// Sampling-free B.2 separation certificate used by the B.19 path.
// A local interval I around the unique central zero of H* is certified to have
// Z*>2 delta by the A.4 pressure/Z enclosure.  Outside I, the explicit cubic
// H* is bounded away from zero using its endpoints and its two critical points.
// sigma is then chosen so chi>chi_threshold outside I.  Consequently
// {chi<=chi_threshold} is contained in I and Z*>2 delta there.
inline AppendixBSeparationAudit appendix_b_separation_interval_audit(
    const OuterProfileParameters& p, double j0 = 0.05,
    double Tf = 20.0, double co = 0.05, double max_step = 0.01,
    double local_radius = 0.02, double chi_threshold = 0.99) {
    p.validate();
    if (!(local_radius > 0.0 && local_radius < 0.2 &&
          chi_threshold > 0.0 && chi_threshold < 1.0))
        throw std::invalid_argument("invalid Appendix B separation parameters");

    const double h = std::exp(p.log_h);
    const double eta0 = appendix_b_find_H_zero(p, j0);
    const double a = std::max(-1.0, eta0 - local_radius);
    const double b = std::min(1.0, eta0 + local_radius);

    // A strictly positive placeholder delta lets the Z enclosure run; the
    // certified delta is selected afterwards from its positive lower bound.
    const auto z = appendix_b_Z_interval_audit(
        p, a, b, std::numeric_limits<double>::min(), j0, Tf, co, max_step);
    if (!z.pressure_enclosure_certified || !z.positive_on_interval ||
        !(z.Z_lower > 0.0))
        throw std::runtime_error("failed to certify Z*>0 near the H* zero");
    const double delta = appendix_b_down(0.5 * z.Z_lower);
    if (!(delta > 0.0)) throw std::runtime_error("nonpositive certified delta");

    // H'=0 has two explicit roots.  On each resulting monotone piece, minima
    // of |H| away from the central root occur at an endpoint.  Since I contains
    // eta0, checking a,b, +/-1 and critical points outside I is sufficient.
    const double disc = 4.0*j0*j0 + 48.0*(4.5-h);
    const double rm = (-2.0*j0-std::sqrt(disc))/24.0;
    const double rp = (-2.0*j0+std::sqrt(disc))/24.0;
    double Hmin = std::min(std::abs(appendix_b_Hstar(a,h,j0)),
                           std::abs(appendix_b_Hstar(b,h,j0)));
    for (double x : {-1.0, rm, rp, 1.0}) {
        if (x < a || x > b)
            Hmin = std::min(Hmin, std::abs(appendix_b_Hstar(x,h,j0)));
    }
    if (!(Hmin > 0.0 && std::isfinite(Hmin)))
        throw std::runtime_error("failed to separate H* outside axial interval");

    // chi=H^2/(H^2+sigma^2)>q iff sigma<|H|sqrt((1-q)/q).
    // Factor 1/2 leaves a strict outward-rounding margin.
    const double sigma = appendix_b_down(
        0.5 * Hmin * std::sqrt((1.0-chi_threshold)/chi_threshold));
    const double chi_min = appendix_b_down(
        Hmin*Hmin/(Hmin*Hmin + sigma*sigma));

    AppendixBSeparationAudit out;
    out.eta0 = eta0;
    out.normalized_Zstar_eta0 = z.Z_center;
    out.normalized_delta_star = delta;
    out.sigma_star = sigma;
    out.min_chi_on_small_Z = chi_min;
    out.unique_H_zero = eta0 > -1.0 && eta0 < 0.0;
    out.positive_Z_at_H_zero = z.Z_lower > 0.0;
    out.chi_separation_certified = chi_min > chi_threshold;
    return out;
}

} // namespace nsblowup
