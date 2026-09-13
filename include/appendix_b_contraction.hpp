#pragma once

#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct AppendixBContractionAudit {
    double inverse_bound{};
    double source_norm{};
    double quadratic_bound{};
    double Lambda{};
    double radius{};
    double image_bound{};
    double lipschitz_bound{};
    bool invariant_ball{};
    bool contraction{};
};

// Banach fixed-point audit for a coefficient-space map of the form
//   z -> L^{-1}(g + Lambda^{-1} Q(z,z)),
// with ||L^{-1}|| <= beta and ||Q(x,y)|| <= q ||x|| ||y||.
//
// This is a reusable bound engine only.  Appendix B is marked complete only
// after the manuscript's actual nonlinear map has been reduced to these
// quantities without replacing any coefficient by an inferred formula.
inline AppendixBContractionAudit appendix_b_quadratic_contraction_audit(
    double beta,
    double source_norm,
    double quadratic_bound,
    double Lambda,
    double radius) {
    if (!(beta >= 0.0) || !(source_norm >= 0.0) ||
        !(quadratic_bound >= 0.0) || !(Lambda > 0.0) || !(radius > 0.0))
        throw std::invalid_argument("invalid Appendix B contraction parameters");

    AppendixBContractionAudit out;
    out.inverse_bound = beta;
    out.source_norm = source_norm;
    out.quadratic_bound = quadratic_bound;
    out.Lambda = Lambda;
    out.radius = radius;
    out.image_bound = beta * (source_norm + quadratic_bound * radius * radius / Lambda);
    out.lipschitz_bound = 2.0 * beta * quadratic_bound * radius / Lambda;
    out.invariant_ball = out.image_bound <= radius;
    out.contraction = out.lipschitz_bound < 1.0;
    return out;
}

// A convenient radius tied to the linear response.  For radius=2 beta ||g||,
// invariance and contraction are both implied by
//   4 beta^2 q ||g|| / Lambda <= 1.
inline AppendixBContractionAudit appendix_b_doubled_linear_ball_audit(
    double beta,
    double source_norm,
    double quadratic_bound,
    double Lambda) {
    const double radius = 2.0 * beta * source_norm;
    if (!(radius > 0.0))
        throw std::invalid_argument("doubled linear ball requires positive source norm and inverse bound");
    return appendix_b_quadratic_contraction_audit(
        beta, source_norm, quadratic_bound, Lambda, radius);
}

} // namespace nsblowup
