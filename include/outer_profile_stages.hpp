#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nsblowup {

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

inline double appendix_a_smooth_step(double y) {
    if (y <= 0.0) return 0.0;
    if (y >= 1.0) return 1.0;
    const double a = std::exp(-1.0 / (y * y));
    const double b = std::exp(-1.0 / ((1.0 - y) * (1.0 - y)));
    return a / (a + b);
}

inline double appendix_a_smooth_step_derivative(double y) {
    if (y <= 0.0 || y >= 1.0) return 0.0;
    const double s = appendix_a_smooth_step(y);
    return s * (1.0 - s) * (2.0 / (y*y*y) + 2.0 / std::pow(1.0-y,3));
}

inline double outer_shape_f(double eta) { return 1.0 / (1.0 + eta * eta); }
inline double outer_reference_U(double eta) { return 4.0 * eta; }
inline double log_outer_reference_E(double x, double eta, double log_Pstar) {
    if (!(x > 0.0 && x <= 1.0)) throw std::invalid_argument("reference x must be in (0,1]");
    return log_Pstar + std::log(outer_shape_f(eta)) + 0.1 * std::log(x);
}
inline double outer_reference_E(double x, double eta, double Pstar) {
    if (!(Pstar > 0.0)) throw std::invalid_argument("Pstar must be positive");
    return std::exp(log_outer_reference_E(x, eta, std::log(Pstar)));
}

inline double outer_l_first_transition(double y) {
    return 0.6 * (1.0 - appendix_a_smooth_step(y));
}
inline double outer_axial_k(double y, double Md) {
    if (y < 0.0) throw std::invalid_argument("axial-stage y must be non-negative");
    return 4.0 * (1.0 - appendix_a_smooth_step(std::log1p(y) / Md));
}
inline double outer_axial_U(double y, double eta, double Md) {
    return outer_axial_k(y, Md) * eta;
}
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
    return {Tw-25.0,Tw-20.0,Tw-20.0,Tw-15.0,Tw-14.0,Tw-9.0,Tw-8.0,Tw-3.0};
}

// Equation (A.10): remove eta dependence of E over a fixed logarithmic interval Tf.
inline double outer_theta_f(double y, double Tf) {
    if (!(Tf > 0.0)) throw std::invalid_argument("Tf must be positive");
    return 1.0 - appendix_a_smooth_step(y / Tf);
}
inline double log_outer_interpolated_E(double y, double eta, double log_eend,
                                       double lambda, double Tf) {
    const double theta = outer_theta_f(y,Tf);
    const double J0 = std::log1p(eta*eta);
    return log_eend - (0.5 + lambda)*y - theta*J0 - (1.0-theta)*std::log(2.0);
}

// Equation (A.12), terminal radial interval 0<=y<=3.
inline double outer_terminal_psi(double y) {
    return 1.0 - appendix_a_smooth_step((y-1.0)/2.0);
}
inline double outer_terminal_f(double y, double rho_o) {
    return 1.0 - rho_o * outer_terminal_psi(y);
}
inline double outer_terminal_fprime(double y, double rho_o) {
    return 0.5 * rho_o * appendix_a_smooth_step_derivative((y-1.0)/2.0);
}
inline double outer_terminal_l(double y, double h, double rho_o) {
    const double f = outer_terminal_f(y,rho_o);
    if (!(f > 0.0)) throw std::invalid_argument("terminal factor must stay positive");
    return -h + outer_terminal_fprime(y,rho_o)/f;
}

// Equation (A.13), evaluated by composite Simpson quadrature.
inline double outer_terminal_Qp(double h, double rho_o, int panels = 1200) {
    if (panels < 2) panels = 2;
    if (panels % 2) ++panels;
    const double dy = 3.0 / panels;
    auto inner_exp = [&](double v) {
        const int n = 200;
        const double ds = v / n;
        double acc = 0.0;
        for (int i=0;i<=n;++i) {
            const double s=i*ds;
            const double w=(i==0||i==n)?1.0:(i%2?4.0:2.0);
            acc += w*(1.0 + outer_terminal_l(s,h,rho_o));
        }
        return std::exp(acc*ds/3.0);
    };
    double sum=0.0;
    for (int i=0;i<=panels;++i) {
        const double v=i*dy;
        const double f=outer_terminal_f(v,rho_o);
        const double integrand=inner_exp(v)*outer_terminal_fprime(v,rho_o)/f;
        const double w=(i==0||i==panels)?1.0:(i%2?4.0:2.0);
        sum += w*integrand;
    }
    return sum*dy/3.0;
}

}  // namespace nsblowup
