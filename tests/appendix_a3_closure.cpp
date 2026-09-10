#include "appendix_a3_closure.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double rel = 1e-8, double abs = 1e-11) {
    return std::abs(a - b) <= abs + rel * std::max({1.0, std::abs(a), std::abs(b)});
}

int fail(const char* msg, double a = 0.0, double b = 0.0) {
    std::cerr << msg << ": " << a << " vs " << b << '\n';
    return 1;
}

}  // namespace

int main() {
    using namespace nsblowup;

    // Pulse shape from A.2, used quantitatively in A.19.
    if (!close(appendix_a3_phi_b(0.0), 0.0)) return fail("phi_b(0)");
    if (!close(appendix_a3_phi_b(0.02), 0.01, 2e-8))
        return fail("phi_b(.02)", appendix_a3_phi_b(0.02), 0.01);
    if (!(appendix_a3_R0(1.0) > 0.0)) return fail("R0 must be positive on main pulse");
    if (!close(appendix_a3_R0(11.0), 0.0)) return fail("R0 must end at xi=11");

    const double Kb = appendix_a3_Kb();
    if (!(Kb > 0.20 && Kb <= 0.25))
        return fail("A.19 Kb bound", Kb, 0.20);

    // Principal A.19 equation has a unique root in the paper's [0.9,1.2] bracket.
    const double amp0 = appendix_a3_solve_amplitude(0.0);
    const double constant = (1.0 - std::exp(-26.0)) / 4.0;
    if (!(amp0 > 0.9 && amp0 < 1.2)) return fail("A.19 amplitude bracket", amp0, 1.0);
    if (!close(amp0 * amp0 * Kb, constant, 2e-9))
        return fail("A.19 principal root residual", amp0 * amp0 * Kb, constant);

    // A small eta-dependent error term keeps the selected root smooth and bracketed.
    auto small_error = [](double amp, double eta) { return 1.0e-4 * eta * amp; };
    const double amp_plus = appendix_a3_solve_amplitude(0.5, small_error);
    const double amp_minus = appendix_a3_solve_amplitude(-0.5, small_error);
    if (!(amp_plus < amp0 && amp_minus > amp0))
        return fail("A.20 smooth amplitude response", amp_plus, amp_minus);

    // A.15: the two late pulse bumps close M and J simultaneously.
    constexpr double lambda = 0.02;
    const auto pulse = appendix_a3_close_pulse_MJ(lambda, 2.0e-7, -3.0e-7, amp0);
    if (!close(pulse.residual_M, 0.0, 0.0, 2e-12))
        return fail("A.15 M closure", pulse.residual_M, 0.0);
    if (!close(pulse.residual_J, 0.0, 0.0, 2e-12))
        return fail("A.15 J closure", pulse.residual_J, 0.0);
    if (!std::isfinite(pulse.c1) || !std::isfinite(pulse.c2))
        return fail("A.15 correction coefficients must be finite");

    // A.11/A.3: relative E bumps prescribe I while preserving pressure increment.
    const auto angular = appendix_a3_close_angular_I_pressure(lambda, 2.0e-4, 0.0);
    if (!close(angular.residual_I, 0.0, 0.0, 2e-11))
        return fail("A.11 I closure", angular.residual_I, 0.0);
    if (!close(angular.residual_pressure, 0.0, 0.0, 2e-11))
        return fail("A.11 pressure closure", angular.residual_pressure, 0.0);
    if (!(std::abs(angular.c1) < 0.1 && std::abs(angular.c2) < 0.1))
        return fail("A.11 small correction branch", angular.c1, angular.c2);

    // A.16 must reproduce Qp at the left endpoint and vanish at the right endpoint.
    constexpr double h = 1.0e-3;
    constexpr double rho_o = 1.0e-4;
    const double Qp = outer_terminal_Qp(h, rho_o, 1000);
    const double Qs0 = appendix_a3_terminal_Qs(0.0, h, rho_o, 1000);
    if (!close(Qs0, Qp, 3e-5, 2e-9))
        return fail("A.16 Qs(0)=Qp", Qs0, Qp);
    if (!close(appendix_a3_terminal_Qs(3.0, h, rho_o), 0.0))
        return fail("A.16 Qs(3)=0");

    // A.17 exact factors for a hold interval of length 4 log(1/h), l=-1.
    const double hold = 4.0 * std::log(1.0 / h);
    const double xe2_factor = std::exp(-2.0 * hold);
    const double e_factor = std::exp(-1.5 * hold);
    if (!close(xe2_factor, std::pow(h, 8), 2e-12))
        return fail("A.17 XE^2 factor", xe2_factor, std::pow(h, 8));
    if (!close(e_factor, std::pow(h, 6), 2e-12))
        return fail("A.17 E factor", e_factor, std::pow(h, 6));

    std::cout << "Appendix A.3 moment closure checks passed: Kb=" << Kb
              << " amplitude=" << amp0 << '\n';
    return 0;
}
