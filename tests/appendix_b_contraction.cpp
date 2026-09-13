#include "appendix_b_contraction.hpp"

#include <cmath>
#include <iostream>

int main() {
    using nsblowup::appendix_b_doubled_linear_ball_audit;
    using nsblowup::appendix_b_quadratic_contraction_audit;

    {
        const auto a = appendix_b_quadratic_contraction_audit(
            2.0, 0.1, 3.0, 100.0, 0.5);
        const double image_expected = 2.0 * (0.1 + 3.0 * 0.25 / 100.0);
        const double lip_expected = 2.0 * 2.0 * 3.0 * 0.5 / 100.0;
        if (!(std::abs(a.image_bound - image_expected) < 1e-15 &&
              std::abs(a.lipschitz_bound - lip_expected) < 1e-15 &&
              a.invariant_ball && a.contraction)) {
            std::cerr << "Appendix B generic contraction audit failed\n";
            return 1;
        }
    }

    {
        const double beta = 1.5;
        const double g = 0.2;
        const double q = 4.0;
        const double critical = 4.0 * beta * beta * q * g;

        const auto good = appendix_b_doubled_linear_ball_audit(
            beta, g, q, 2.0 * critical);
        if (!(good.invariant_ball && good.contraction && good.lipschitz_bound < 1.0)) {
            std::cerr << "Appendix B doubled-linear ball should contract above threshold\n";
            return 2;
        }

        const auto bad = appendix_b_doubled_linear_ball_audit(
            beta, g, q, 0.5 * critical);
        if (bad.invariant_ball || bad.contraction) {
            std::cerr << "Appendix B doubled-linear ball threshold regression failed\n";
            return 3;
        }
    }

    return 0;
}
