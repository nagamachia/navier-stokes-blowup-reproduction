#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

namespace {
using Complex = std::complex<double>;
using Vec3 = std::array<Complex, 3>;

int wave_number(int index, int n) {
    return (index <= n / 2) ? index : index - n;
}

Vec3 project_mode(const Vec3& velocity, const std::array<int, 3>& k) {
    const double k2 = static_cast<double>(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
    if (k2 == 0.0) {
        return velocity;
    }

    const Complex k_dot_u = static_cast<double>(k[0]) * velocity[0]
                          + static_cast<double>(k[1]) * velocity[1]
                          + static_cast<double>(k[2]) * velocity[2];

    Vec3 projected = velocity;
    for (int component = 0; component < 3; ++component) {
        projected[component] -= static_cast<double>(k[component]) * k_dot_u / k2;
    }
    return projected;
}

double norm_squared(const Vec3& value) {
    return std::norm(value[0]) + std::norm(value[1]) + std::norm(value[2]);
}

double max_component_difference(const Vec3& a, const Vec3& b) {
    return std::max({std::abs(a[0] - b[0]), std::abs(a[1] - b[1]), std::abs(a[2] - b[2])});
}
}  // namespace

int main() {
    constexpr int n = 8;
    constexpr double tolerance = 1e-12;

    double max_divergence = 0.0;
    double max_idempotence_error = 0.0;
    double max_energy_increase = 0.0;
    std::size_t tested_modes = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int l = 0; l < n; ++l) {
                const std::array<int, 3> k{
                    wave_number(i, n),
                    wave_number(j, n),
                    wave_number(l, n)};

                if (k[0] == 0 && k[1] == 0 && k[2] == 0) {
                    continue;
                }

                // Deterministic complex coefficients exercise both real and imaginary parts.
                const Vec3 velocity{
                    Complex{0.25 + 0.11 * i, -0.17 + 0.07 * j},
                    Complex{-0.31 + 0.05 * j, 0.23 - 0.03 * l},
                    Complex{0.19 - 0.04 * l, -0.29 + 0.02 * i}};

                const Vec3 projected = project_mode(velocity, k);
                const Complex divergence = static_cast<double>(k[0]) * projected[0]
                                         + static_cast<double>(k[1]) * projected[1]
                                         + static_cast<double>(k[2]) * projected[2];
                max_divergence = std::max(max_divergence, std::abs(divergence));

                const Vec3 projected_twice = project_mode(projected, k);
                max_idempotence_error = std::max(
                    max_idempotence_error,
                    max_component_difference(projected, projected_twice));

                const double energy_increase = norm_squared(projected) - norm_squared(velocity);
                max_energy_increase = std::max(max_energy_increase, energy_increase);
                ++tested_modes;
            }
        }
    }

    std::cout << std::setprecision(16)
              << "N=" << n << '\n'
              << "tested_modes=" << tested_modes << '\n'
              << "max_abs_k_dot_u_after_projection=" << max_divergence << '\n'
              << "max_idempotence_error=" << max_idempotence_error << '\n'
              << "max_energy_increase=" << max_energy_increase << '\n';

    if (max_divergence > tolerance ||
        max_idempotence_error > tolerance ||
        max_energy_increase > tolerance) {
        std::cerr << "Divergence-free projection check failed\n";
        return 1;
    }

    std::cout << "Divergence-free projection check passed\n";
    return 0;
}
