#pragma once

#include "leading_profile_relations.hpp"

#include <cmath>
#include <stdexcept>

namespace nsblowup {

struct LeadingStressPoint {
    double H;
    double F;
    double ell;
    double W;
    double Hc;
    double S_q;
    double S_n;
    double Q_s;
    double N_s;
    double p_s_theta;
    double p_s_z;
    double a;
    double b_s;
    double T0_theta;
    double T0_z;
};

// Profile callbacks use the similarity variables (X, eta):
//   F, D_X F, partial_eta F,
//   U, D_X U, partial_eta U,
//   Pi, D_X Pi, partial_eta Pi.
// E is reconstructed from the regular-axis factorization E=sqrt(2X)F.
template <
    class FFunc,
    class DXFFunc,
    class DEtaFFunc,
    class UFunc,
    class DXUFunc,
    class DEtaUFunc,
    class PiFunc,
    class DXPiFunc,
    class DEtaPiFunc>
LeadingStressPoint leading_stress_point(
    FFunc&& F_fn,
    DXFFunc&& DXF_fn,
    DEtaFFunc&& dF_deta_fn,
    UFunc&& U_fn,
    DXUFunc&& DXU_fn,
    DEtaUFunc&& dU_deta_fn,
    PiFunc&& Pi_fn,
    DXPiFunc&& DXPi_fn,
    DEtaPiFunc&& dPi_deta_fn,
    double X,
    double eta,
    double h,
    int panels = 256) {
    if (X < 0.0) {
        throw std::invalid_argument("X must be non-negative");
    }
    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }

    const double A = 0.5 + h;
    const double D = 0.5 - h;
    const double d = 1.0 - eta * eta;
    const double L = 1.0 - 2.0 * h * eta * eta;

    const auto average_U = [&](double x, double e) {
        return radial_average(U_fn, x, e, panels);
    };
    const auto d_eta_average_U = [&](double x, double e) {
        return radial_average(dU_deta_fn, x, e, panels);
    };

    auto source_q = [&](double x, double e) {
        const double F = F_fn(x, e);
        if (!(F > 0.0)) {
            throw std::runtime_error("F must be positive for log derivatives");
        }
        const double ell = 1.0 + DXF_fn(x, e) / F;
        const double dd = 1.0 - e * e;
        const double LL = 1.0 - 2.0 * h * e * e;
        (void)LL;
        const double W = 1.0
            - 2.0 * D * e * average_U(x, e)
            - dd * d_eta_average_U(x, e);
        const double Hc = D * e + dd * U_fn(x, e);
        const double logE_eta = dF_deta_fn(x, e) / F;
        return -W * ell
            - h * (1.0 - 2.0 * e * U_fn(x, e))
            - Hc * logE_eta;
    };

    auto source_n = [&](double x, double e) {
        const double dd = 1.0 - e * e;
        const double W = 1.0
            - 2.0 * D * e * average_U(x, e)
            - dd * d_eta_average_U(x, e);
        const double Hc = D * e + dd * U_fn(x, e);
        return -W * DXU_fn(x, e)
            - A * (1.0 - 2.0 * e * U_fn(x, e)) * U_fn(x, e)
            - Hc * dU_deta_fn(x, e)
            - dd * dPi_deta_fn(x, e)
            + 4.0 * A * e * Pi_fn(x, e)
            + 2.0 * e * DXPi_fn(x, e);
    };

    const double F = F_fn(X, eta);
    if (!(F > 0.0)) {
        throw std::runtime_error("F must be positive at evaluation point");
    }

    const double DXF = DXF_fn(X, eta);
    const double ell = 1.0 + DXF / F;
    const double averageU = average_U(X, eta);
    const double dEtaAverageU = d_eta_average_U(X, eta);
    const double W = 1.0 - 2.0 * D * eta * averageU - d * dEtaAverageU;
    const double Hc = D * eta + d * U_fn(X, eta);
    const double Sq = source_q(X, eta);
    const double Sn = source_n(X, eta);

    // Equation (4.10), written in the axis-regular s in [0,1] form.
    const double Qs = simpson_integral(
        [&](double s) {
            return s * F_fn(s * X, eta) * source_q(s * X, eta);
        },
        0.0,
        1.0,
        panels) / F;

    const double Ns = simpson_integral(
        [&](double s) { return source_n(s * X, eta); },
        0.0,
        1.0,
        panels);

    const double E = X == 0.0 ? 0.0 : std::sqrt(2.0 * X) * F;
    const double H = 2.0 * X * F;
    const double a = 2.0 - 2.0 * ell;

    double bs = 0.0;
    double ps_z = 0.0;
    if (X > 0.0) {
        bs = 2.0 * DXU_fn(X, eta) / E;
        ps_z = X * Ns / (L * E);
    }
    const double ps_theta = X * Qs / L;

    // Equation (4.11): T0 = F (p_s - (a,-b_s)).
    const double T0_theta = F * (ps_theta - a);
    const double T0_z = F * (ps_z + bs);

    return {
        H,
        F,
        ell,
        W,
        Hc,
        Sq,
        Sn,
        Qs,
        Ns,
        ps_theta,
        ps_z,
        a,
        bs,
        T0_theta,
        T0_z,
    };
}

}  // namespace nsblowup
