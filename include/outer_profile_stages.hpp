#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

// Appendix A.2 has the hierarchy Td=exp(Md)+10, P*>exp(Td),
// h<exp(-Td). P* and h therefore cannot in general be stored directly
// in IEEE double. Keep their logarithms for schedule validation.
struct OuterProfileParameters {
    double Md{4.0};
    double log_Pstar{70.0};
    double lambda{0.02};
    double log_h{-80.0};

    double Td() const { return std::exp(Md) + 10.0; }
    double Tw() const { return 60.0 * std::log(1.0 / lambda); }

    void validate() const {
        if (!(Md > 0.0)) throw std::invalid_argument("Md must be positive");
        if (!(lambda > 0.0 && lambda < 1.0)) throw std::invalid_argument("lambda must lie in (0,1)");
        const double td = Td();
        if (!(log_Pstar > td)) throw std::invalid_argument("A.6 requires log(Pstar) > Td");
        if (!(log_h < std::min({std::log(0.01), std::log(lambda), -td})))
            throw std::invalid_argument("A.6 requires log(h) < min(log(1/100),log(lambda),-Td)");
    }
};

// Equation (A.5).
inline double appendix_a_smooth_step(double y) {
    if (y <= 0.0) return 0.0;
    if (y >= 1.0) return 1.0;
    const double a = std::exp(-1.0 / (y * y));
    const double b = std::exp(-1.0 / ((1.0 - y) * (1.0 - y)));
    return a / (a + b);
}

inline double outer_shape_f(double eta) { return 1.0 / (1.0 + eta * eta); }

// Equation (A.7), x=X/XR <= 1. A logarithmic form is supplied for
// theorem-scale amplitudes; the direct form is convenient for normalized tests.
inline double outer_reference_U(double eta) { return 4.0 * eta; }
inline double log_outer_reference_E(double x, double eta, double log_Pstar) {
    if (!(x > 0.0 && x <= 1.0)) throw std::invalid_argument("reference x must be in (0,1]");
    return log_Pstar + std::log(outer_shape_f(eta)) + 0.1 * std::log(x);
}
inline double outer_reference_E(double x, double eta, double Pstar) {
    if (!(Pstar > 0.0)) throw std::invalid_argument("Pstar must be positive");
    return std::exp(log_outer_reference_E(x, eta, std::log(Pstar)));
}

// First logarithmic stage after x=1: l=(3/5)(1-sigma(y)), 0<=y<=1.
inline double outer_l_first_transition(double y) {
    return 0.6 * (1.0 - appendix_a_smooth_step(y));
}

// Axial-decay stage from A.2: l=0 and U=k(y) eta, 0<=y<=Td.
inline double outer_axial_k(double y, double Md) {
    if (y < 0.0) throw std::invalid_argument("axial-stage y must be non-negative");
    return 4.0 * (1.0 - appendix_a_smooth_step(std::log1p(y) / Md));
}
inline double outer_axial_U(double y, double eta, double Md) {
    return outer_axial_k(y, Md) * eta;
}

// Intermediate power-law stage, equation (A.9).
inline double outer_l_intermediate_entry(double y, double lambda) {
    return -lambda * appendix_a_smooth_step(y);
}
inline double outer_l_intermediate(double lambda) { return -lambda; }

struct ReservedPatchesA9 {
    double profile_a, profile_b;
    double heat_a, heat_b;
    double positive_a, positive_b;
    double mean_a, mean_b;
};

inline ReservedPatchesA9 reserved_patches_a9(double lambda) {
    if (!(lambda > 0.0 && lambda < 1.0)) throw std::invalid_argument("lambda must lie in (0,1)");
    const double Tw = 60.0 * std::log(1.0 / lambda);
    if (!(Tw > 25.0)) throw std::invalid_argument("lambda too large for A.9 reserved patches");
    return {Tw - 25.0, Tw - 20.0,
            Tw - 20.0, Tw - 15.0,
            Tw - 14.0, Tw - 9.0,
            Tw - 8.0, Tw - 3.0};
}

}  // namespace nsblowup
