#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

namespace {
using Complex = std::complex<double>;
using RealField = std::vector<double>;
using SpectralField = std::vector<Complex>;
using RealVectorField = std::array<RealField, 3>;
using SpectralVectorField = std::array<SpectralField, 3>;

struct Grid {
    int n;
    int size() const { return n * n * n; }
    int index(int i, int j, int k) const { return (i * n + j) * n + k; }
    int wave_number(int index_value) const {
        return (index_value <= n / 2) ? index_value : index_value - n;
    }
};

SpectralField dft3(const Grid& grid, const RealField& values) {
    const double two_pi = 2.0 * std::numbers::pi;
    SpectralField spectrum(grid.size(), {0.0, 0.0});
    for (int a = 0; a < grid.n; ++a) {
        for (int b = 0; b < grid.n; ++b) {
            for (int c = 0; c < grid.n; ++c) {
                Complex sum{0.0, 0.0};
                for (int i = 0; i < grid.n; ++i) {
                    for (int j = 0; j < grid.n; ++j) {
                        for (int k = 0; k < grid.n; ++k) {
                            const double phase = -two_pi *
                                static_cast<double>(a * i + b * j + c * k) / grid.n;
                            sum += values[grid.index(i, j, k)] *
                                   Complex{std::cos(phase), std::sin(phase)};
                        }
                    }
                }
                spectrum[grid.index(a, b, c)] =
                    sum / static_cast<double>(grid.size());
            }
        }
    }
    return spectrum;
}

RealField idft3(const Grid& grid, const SpectralField& spectrum) {
    const double two_pi = 2.0 * std::numbers::pi;
    RealField values(grid.size(), 0.0);
    for (int i = 0; i < grid.n; ++i) {
        for (int j = 0; j < grid.n; ++j) {
            for (int k = 0; k < grid.n; ++k) {
                Complex sum{0.0, 0.0};
                for (int a = 0; a < grid.n; ++a) {
                    for (int b = 0; b < grid.n; ++b) {
                        for (int c = 0; c < grid.n; ++c) {
                            const double phase = two_pi *
                                static_cast<double>(a * i + b * j + c * k) / grid.n;
                            sum += spectrum[grid.index(a, b, c)] *
                                   Complex{std::cos(phase), std::sin(phase)};
                        }
                    }
                }
                values[grid.index(i, j, k)] = sum.real();
            }
        }
    }
    return values;
}

void project_and_dealias(const Grid& grid, SpectralVectorField& field) {
    const int cutoff = grid.n / 3;
    for (int a = 0; a < grid.n; ++a) {
        const int kx = grid.wave_number(a);
        for (int b = 0; b < grid.n; ++b) {
            const int ky = grid.wave_number(b);
            for (int c = 0; c < grid.n; ++c) {
                const int kz = grid.wave_number(c);
                const int idx = grid.index(a, b, c);
                if (std::abs(kx) > cutoff ||
                    std::abs(ky) > cutoff ||
                    std::abs(kz) > cutoff) {
                    for (auto& component : field) {
                        component[idx] = {0.0, 0.0};
                    }
                    continue;
                }

                const double k2 =
                    static_cast<double>(kx * kx + ky * ky + kz * kz);
                if (k2 == 0.0) {
                    continue;
                }

                const Complex dot = static_cast<double>(kx) * field[0][idx]
                                  + static_cast<double>(ky) * field[1][idx]
                                  + static_cast<double>(kz) * field[2][idx];
                field[0][idx] -= static_cast<double>(kx) * dot / k2;
                field[1][idx] -= static_cast<double>(ky) * dot / k2;
                field[2][idx] -= static_cast<double>(kz) * dot / k2;
            }
        }
    }
}

SpectralVectorField rhs(
    const Grid& grid,
    const SpectralVectorField& uhat,
    double viscosity) {
    SpectralVectorField omega_hat;
    SpectralVectorField nonlinear_hat;
    for (auto& field : omega_hat) {
        field.assign(grid.size(), {0.0, 0.0});
    }
    for (auto& field : nonlinear_hat) {
        field.assign(grid.size(), {0.0, 0.0});
    }

    const Complex imaginary_unit{0.0, 1.0};
    for (int a = 0; a < grid.n; ++a) {
        const double kx = grid.wave_number(a);
        for (int b = 0; b < grid.n; ++b) {
            const double ky = grid.wave_number(b);
            for (int c = 0; c < grid.n; ++c) {
                const double kz = grid.wave_number(c);
                const int idx = grid.index(a, b, c);
                omega_hat[0][idx] = imaginary_unit *
                    (ky * uhat[2][idx] - kz * uhat[1][idx]);
                omega_hat[1][idx] = imaginary_unit *
                    (kz * uhat[0][idx] - kx * uhat[2][idx]);
                omega_hat[2][idx] = imaginary_unit *
                    (kx * uhat[1][idx] - ky * uhat[0][idx]);
            }
        }
    }

    RealVectorField velocity;
    RealVectorField vorticity;
    for (int component = 0; component < 3; ++component) {
        velocity[component] = idft3(grid, uhat[component]);
        vorticity[component] = idft3(grid, omega_hat[component]);
    }

    RealVectorField cross;
    for (auto& field : cross) {
        field.assign(grid.size(), 0.0);
    }
    for (int idx = 0; idx < grid.size(); ++idx) {
        cross[0][idx] = velocity[1][idx] * vorticity[2][idx]
                      - velocity[2][idx] * vorticity[1][idx];
        cross[1][idx] = velocity[2][idx] * vorticity[0][idx]
                      - velocity[0][idx] * vorticity[2][idx];
        cross[2][idx] = velocity[0][idx] * vorticity[1][idx]
                      - velocity[1][idx] * vorticity[0][idx];
    }

    for (int component = 0; component < 3; ++component) {
        nonlinear_hat[component] = dft3(grid, cross[component]);
    }
    project_and_dealias(grid, nonlinear_hat);

    for (int a = 0; a < grid.n; ++a) {
        const int kx = grid.wave_number(a);
        for (int b = 0; b < grid.n; ++b) {
            const int ky = grid.wave_number(b);
            for (int c = 0; c < grid.n; ++c) {
                const int kz = grid.wave_number(c);
                const int idx = grid.index(a, b, c);
                const double k2 =
                    static_cast<double>(kx * kx + ky * ky + kz * kz);
                for (int component = 0; component < 3; ++component) {
                    nonlinear_hat[component][idx] -=
                        viscosity * k2 * uhat[component][idx];
                }
            }
        }
    }
    return nonlinear_hat;
}

SpectralVectorField add_scaled(
    const SpectralVectorField& a,
    const SpectralVectorField& b,
    double scale) {
    SpectralVectorField out = a;
    for (int component = 0; component < 3; ++component) {
        for (std::size_t i = 0; i < out[component].size(); ++i) {
            out[component][i] += scale * b[component][i];
        }
    }
    return out;
}

SpectralVectorField rk4_step(
    const Grid& grid,
    const SpectralVectorField& velocity,
    double dt,
    double viscosity) {
    const auto k1 = rhs(grid, velocity, viscosity);
    const auto k2 = rhs(grid, add_scaled(velocity, k1, 0.5 * dt), viscosity);
    const auto k3 = rhs(grid, add_scaled(velocity, k2, 0.5 * dt), viscosity);
    const auto k4 = rhs(grid, add_scaled(velocity, k3, dt), viscosity);

    SpectralVectorField out = velocity;
    for (int component = 0; component < 3; ++component) {
        for (std::size_t i = 0; i < out[component].size(); ++i) {
            out[component][i] += (dt / 6.0) *
                (k1[component][i]
                 + 2.0 * k2[component][i]
                 + 2.0 * k3[component][i]
                 + k4[component][i]);
        }
    }
    project_and_dealias(grid, out);
    return out;
}

double max_divergence(
    const Grid& grid,
    const SpectralVectorField& field) {
    double result = 0.0;
    for (int a = 0; a < grid.n; ++a) {
        const int kx = grid.wave_number(a);
        for (int b = 0; b < grid.n; ++b) {
            const int ky = grid.wave_number(b);
            for (int c = 0; c < grid.n; ++c) {
                const int kz = grid.wave_number(c);
                const int idx = grid.index(a, b, c);
                const Complex divergence =
                    static_cast<double>(kx) * field[0][idx]
                    + static_cast<double>(ky) * field[1][idx]
                    + static_cast<double>(kz) * field[2][idx];
                result = std::max(result, std::abs(divergence));
            }
        }
    }
    return result;
}
}  // namespace

int main() {
    const Grid grid{6};
    constexpr double viscosity = 0.1;
    constexpr double dt = 0.01;
    constexpr int steps = 10;
    constexpr double tolerance = 2e-10;
    const double two_pi = 2.0 * std::numbers::pi;

    RealVectorField initial;
    for (auto& field : initial) {
        field.assign(grid.size(), 0.0);
    }

    // Two-dimensional Taylor-Green vortex embedded in the 3D periodic box.
    // Its projected nonlinear term vanishes, so the exact velocity decays as
    // exp(-2 * nu * t), making it a compact end-to-end solver regression.
    for (int i = 0; i < grid.n; ++i) {
        const double x = two_pi * static_cast<double>(i) / grid.n;
        for (int j = 0; j < grid.n; ++j) {
            const double y = two_pi * static_cast<double>(j) / grid.n;
            for (int k = 0; k < grid.n; ++k) {
                const int idx = grid.index(i, j, k);
                initial[0][idx] = std::sin(x) * std::cos(y);
                initial[1][idx] = -std::cos(x) * std::sin(y);
            }
        }
    }

    SpectralVectorField velocity_hat;
    for (int component = 0; component < 3; ++component) {
        velocity_hat[component] = dft3(grid, initial[component]);
    }
    project_and_dealias(grid, velocity_hat);

    for (int step = 0; step < steps; ++step) {
        velocity_hat = rk4_step(grid, velocity_hat, dt, viscosity);
    }

    RealVectorField numerical;
    for (int component = 0; component < 3; ++component) {
        numerical[component] = idft3(grid, velocity_hat[component]);
    }

    const double time = steps * dt;
    const double exact_factor = std::exp(-2.0 * viscosity * time);
    double max_velocity_error = 0.0;
    double initial_energy = 0.0;
    double final_energy = 0.0;
    for (int component = 0; component < 3; ++component) {
        for (int idx = 0; idx < grid.size(); ++idx) {
            const double exact = exact_factor * initial[component][idx];
            max_velocity_error = std::max(
                max_velocity_error,
                std::abs(numerical[component][idx] - exact));
            initial_energy += initial[component][idx] * initial[component][idx];
            final_energy += numerical[component][idx] * numerical[component][idx];
        }
    }

    const double energy_ratio = final_energy / initial_energy;
    const double exact_energy_ratio = std::exp(-4.0 * viscosity * time);
    const double energy_error = std::abs(energy_ratio - exact_energy_ratio);
    const double divergence_error = max_divergence(grid, velocity_hat);

    std::cout << std::setprecision(16)
              << "N=" << grid.n << '\n'
              << "time=" << time << '\n'
              << "max_velocity_error=" << max_velocity_error << '\n'
              << "divergence_error=" << divergence_error << '\n'
              << "energy_ratio=" << energy_ratio << '\n'
              << "exact_energy_ratio=" << exact_energy_ratio << '\n'
              << "energy_error=" << energy_error << '\n';

    if (max_velocity_error > tolerance ||
        divergence_error > tolerance ||
        energy_error > tolerance) {
        std::cerr << "Taylor-Green regression failed\n";
        return 1;
    }

    std::cout << "Taylor-Green regression passed\n";
    return 0;
}
