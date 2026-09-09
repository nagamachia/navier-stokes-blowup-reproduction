#include "similarity_coordinates.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double rel = 2e-13) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel * scale;
}

int fail(const char* message) {
    std::cerr << message << '\n';
    return 1;
}

}  // namespace

int main() {
    constexpr double h = 0.005;
    constexpr double D = 0.5 - h;

    for (double tau : {1.0, 1e-2, 1e-4, 1e-6}) {
        for (double eta : {-0.9, -0.5, 0.0, 0.4, 0.9}) {
            const double q_reference = tau / (1.0 - eta * eta);
            const double z = std::pow(q_reference, D) * eta;
            const double r = std::sqrt(2.0 * 0.7 * q_reference);
            const auto c = nsblowup::coordinates(r, z, tau, h);

            if (!close(c.q, q_reference)) return fail("q identity failed");
            if (!close(c.eta, eta)) return fail("eta identity failed");
            if (!close(c.X, 0.7)) return fail("X identity failed");

            const double reconstructed_tau = c.q * (1.0 - c.eta * c.eta);
            const double reconstructed_z = std::pow(c.q, D) * c.eta;
            if (!close(reconstructed_tau, tau)) return fail("tau reconstruction failed");
            if (!close(reconstructed_z, z)) return fail("z reconstruction failed");
        }
    }

    std::cout << "similarity-coordinate regression checks passed\n";
    return 0;
}
