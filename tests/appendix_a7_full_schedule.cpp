#include "appendix_a7_full_schedule.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

double max_abs(const nsblowup::AppendixA7MomentChange& q) {
    return std::max({std::abs(q.Cp), std::abs(q.S), std::abs(q.I)});
}

} // namespace

int main() {
    using namespace nsblowup;

    constexpr double eta = 0.35;
    constexpr double lambda = 0.05;
    constexpr double h = 0.08;
    constexpr double XK = 120.0;

    const auto d = appendix_a7_heat_tail_discrepancy_normalized(
        eta, h, XK, 12.0, 512);
    if (!(std::isfinite(d.Cp) && std::isfinite(d.S) && std::isfinite(d.I))) {
        std::cerr << "non-finite A.7 heat-tail discrepancy\n";
        return 1;
    }
    if (!(d.Cp < 0.0 && d.S > 0.0 && d.I < 0.0)) {
        std::cerr << "unexpected A.7 heat-tail discrepancy signs\n";
        return 2;
    }

    const auto coupled = appendix_a7_couple_heat_to_reserved_patch(
        eta, lambda, h, XK, 1.0, 1.0, 12.0, 512);
    const double scale = std::max(1e-14, max_abs(coupled.heat_discrepancy));
    if (max_abs(coupled.residual) > 2e-9 * scale) {
        std::cerr << "A.7 replacement/compensation residual too large: "
                  << max_abs(coupled.residual) << " scale=" << scale << "\n";
        return 3;
    }

    const auto plus = appendix_a7_scale_change({1.0, -2.0, 3.0}, 4.0, 5.0);
    const auto back = appendix_a7_unscale_change(plus, 4.0, 5.0);
    if (std::abs(back.Cp - 1.0) > 1e-13 ||
        std::abs(back.S + 2.0) > 1e-13 ||
        std::abs(back.I - 3.0) > 1e-13) {
        std::cerr << "A.7 scale conversion is not invertible\n";
        return 4;
    }

    return 0;
}
