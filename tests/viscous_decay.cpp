#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

namespace {
using Complex = std::complex<double>;
using Vec3 = std::array<Complex, 3>;

Vec3 rhs(const Vec3& u, double decay_rate) {
    return Vec3{-decay_rate * u[0], -decay_rate * u[1], -decay_rate * u[2]};
}

Vec3 add_scaled(const Vec3& a, const Vec3& b, double scale) {
    return Vec3{a[0] + scale * b[0], a[1] + scale * b[1], a[2] + scale * b[2]};
}

Vec3 rk4_step(const Vec3& u, double dt, double decay_rate) {
    const Vec3 k1 = rhs(u, decay_rate);
    const Vec3 k2 = rhs(add_scaled(u, k1, 0.5 * dt), decay_rate);
    const Vec3 k3 = rhs(add_scaled(u, k2, 0.5 * dt), decay_rate);
    const Vec3 k4 = rhs(add_scaled(u, k3, dt), decay_rate);

    Vec3 out = u;
    for (int c = 0; c < 3; ++c) {
        out[c] += (dt / 6.0) * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
    }
    return out;
}

double norm_squared(const Vec3& u) {
    return std::norm(u[0]) + std::norm(u[1]) + std::norm(u[2]);
}
}  // namespace

int main() {
    constexpr std::array<int, 3> k{2, 1, 0};
    constexpr double viscosity = 0.1;
    constexpr double final_time = 1.0;
    constexpr double dt = 1e-3;
    constexpr double tolerance = 1e-12;

    const double k2 = static_cast<double>(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
    const double decay_rate = viscosity * k2;

    // k dot u = 2*(1) + 1*(-2) = 0, so this is an incompressible Fourier mode.
    const Vec3 initial{Complex{1.0, 0.25}, Complex{-2.0, -0.5}, Complex{0.3, -0.4}};
    Vec3 numerical = initial;

    const int steps = static_cast<int>(std::lround(final_time / dt));
    for (int step = 0; step < steps; ++step) {
        numerical = rk4_step(numerical, dt, decay_rate);
    }

    const double exact_factor = std::exp(-decay_rate * final_time);
    Vec3 exact = initial;
    for (Complex& value : exact) {
        value *= exact_factor;
    }

    double max_error = 0.0;
    for (int c = 0; c < 3; ++c) {
        max_error = std::max(max_error, std::abs(numerical[c] - exact[c]));
    }

    const double energy_ratio = norm_squared(numerical) / norm_squared(initial);
    const double exact_energy_ratio = std::exp(-2.0 * decay_rate * final_time);
    const double energy_error = std::abs(energy_ratio - exact_energy_ratio);

    std::cout << std::setprecision(16)
              << "wave_number_squared=" << k2 << '\n'
              << "viscosity=" << viscosity << '\n'
              << "steps=" << steps << '\n'
              << "max_mode_error=" << max_error << '\n'
              << "energy_ratio=" << energy_ratio << '\n'
              << "exact_energy_ratio=" << exact_energy_ratio << '\n'
              << "energy_error=" << energy_error << '\n';

    if (max_error > tolerance || energy_error > tolerance) {
        std::cerr << "Viscous Fourier-mode decay check failed\n";
        return 1;
    }

    std::cout << "Viscous Fourier-mode decay check passed\n";
    return 0;
}
