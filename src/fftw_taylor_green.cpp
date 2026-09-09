#include <fftw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;
using RealField = std::vector<double>;
using SpectralField = std::vector<Complex>;
using RealVectorField = std::array<RealField, 3>;
using SpectralVectorField = std::array<SpectralField, 3>;

struct Grid {
    int n;

    std::size_t size() const {
        return static_cast<std::size_t>(n) * n * n;
    }

    std::size_t index(int i, int j, int k) const {
        return (static_cast<std::size_t>(i) * n + j) * n + k;
    }

    int wave_number(int index_value) const {
        return (index_value <= n / 2) ? index_value : index_value - n;
    }
};

class Fft3d {
public:
    explicit Fft3d(const Grid& grid)
        : grid_(grid),
          input_(static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * grid.size()))),
          output_(static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * grid.size()))) {
        if (input_ == nullptr || output_ == nullptr) {
            throw std::runtime_error("failed to allocate FFTW buffers");
        }

        forward_plan_ = fftw_plan_dft_3d(
            grid.n, grid.n, grid.n, input_, output_, FFTW_FORWARD, FFTW_ESTIMATE);
        backward_plan_ = fftw_plan_dft_3d(
            grid.n, grid.n, grid.n, input_, output_, FFTW_BACKWARD, FFTW_ESTIMATE);

        if (forward_plan_ == nullptr || backward_plan_ == nullptr) {
            throw std::runtime_error("failed to create FFTW plans");
        }
    }

    Fft3d(const Fft3d&) = delete;
    Fft3d& operator=(const Fft3d&) = delete;

    ~Fft3d() {
        if (forward_plan_ != nullptr) {
            fftw_destroy_plan(forward_plan_);
        }
        if (backward_plan_ != nullptr) {
            fftw_destroy_plan(backward_plan_);
        }
        fftw_free(input_);
        fftw_free(output_);
    }

    SpectralField forward(const RealField& values) {
        if (values.size() != grid_.size()) {
            throw std::invalid_argument("forward FFT field size mismatch");
        }

        for (std::size_t i = 0; i < grid_.size(); ++i) {
            input_[i][0] = values[i];
            input_[i][1] = 0.0;
        }
        fftw_execute(forward_plan_);

        const double normalization = 1.0 / static_cast<double>(grid_.size());
        SpectralField spectrum(grid_.size());
        for (std::size_t i = 0; i < grid_.size(); ++i) {
            spectrum[i] = normalization * Complex{output_[i][0], output_[i][1]};
        }
        return spectrum;
    }

    RealField inverse(const SpectralField& spectrum) {
        if (spectrum.size() != grid_.size()) {
            throw std::invalid_argument("inverse FFT field size mismatch");
        }

        for (std::size_t i = 0; i < grid_.size(); ++i) {
            input_[i][0] = spectrum[i].real();
            input_[i][1] = spectrum[i].imag();
        }
        fftw_execute(backward_plan_);

        RealField values(grid_.size());
        for (std::size_t i = 0; i < grid_.size(); ++i) {
            values[i] = output_[i][0];
        }
        return values;
    }

private:
    Grid grid_;
    fftw_complex* input_{nullptr};
    fftw_complex* output_{nullptr};
    fftw_plan forward_plan_{nullptr};
    fftw_plan backward_plan_{nullptr};
};

SpectralVectorField make_spectral_vector(std::size_t size) {
    SpectralVectorField field;
    for (auto& component : field) {
        component.assign(size, {0.0, 0.0});
    }
    return field;
}

RealVectorField make_real_vector(std::size_t size) {
    RealVectorField field;
    for (auto& component : field) {
        component.assign(size, 0.0);
    }
    return field;
}

void project_and_dealias(const Grid& grid, SpectralVectorField& field) {
    const int cutoff = grid.n / 3;
    for (int i = 0; i < grid.n; ++i) {
        const int kx = grid.wave_number(i);
        for (int j = 0; j < grid.n; ++j) {
            const int ky = grid.wave_number(j);
            for (int k = 0; k < grid.n; ++k) {
                const int kz = grid.wave_number(k);
                const std::size_t idx = grid.index(i, j, k);

                if (std::abs(kx) > cutoff ||
                    std::abs(ky) > cutoff ||
                    std::abs(kz) > cutoff) {
                    for (auto& component : field) {
                        component[idx] = {0.0, 0.0};
                    }
                    continue;
                }

                const double k_squared =
                    static_cast<double>(kx * kx + ky * ky + kz * kz);
                if (k_squared == 0.0) {
                    continue;
                }

                const Complex k_dot_u =
                    static_cast<double>(kx) * field[0][idx]
                    + static_cast<double>(ky) * field[1][idx]
                    + static_cast<double>(kz) * field[2][idx];

                field[0][idx] -= static_cast<double>(kx) * k_dot_u / k_squared;
                field[1][idx] -= static_cast<double>(ky) * k_dot_u / k_squared;
                field[2][idx] -= static_cast<double>(kz) * k_dot_u / k_squared;
            }
        }
    }
}

SpectralVectorField vorticity_hat(
    const Grid& grid,
    const SpectralVectorField& velocity_hat) {
    SpectralVectorField omega_hat = make_spectral_vector(grid.size());
    const Complex imaginary_unit{0.0, 1.0};

    for (int i = 0; i < grid.n; ++i) {
        const double kx = grid.wave_number(i);
        for (int j = 0; j < grid.n; ++j) {
            const double ky = grid.wave_number(j);
            for (int k = 0; k < grid.n; ++k) {
                const double kz = grid.wave_number(k);
                const std::size_t idx = grid.index(i, j, k);
                omega_hat[0][idx] = imaginary_unit *
                    (ky * velocity_hat[2][idx] - kz * velocity_hat[1][idx]);
                omega_hat[1][idx] = imaginary_unit *
                    (kz * velocity_hat[0][idx] - kx * velocity_hat[2][idx]);
                omega_hat[2][idx] = imaginary_unit *
                    (kx * velocity_hat[1][idx] - ky * velocity_hat[0][idx]);
            }
        }
    }
    return omega_hat;
}

double spectral_l2(const SpectralVectorField& field) {
    double sum = 0.0;
    for (const auto& component : field) {
        for (const Complex& value : component) {
            sum += std::norm(value);
        }
    }
    return std::sqrt(sum);
}

SpectralVectorField rhs(
    const Grid& grid,
    Fft3d& fft,
    const SpectralVectorField& velocity_hat,
    double viscosity,
    double* projected_nonlinear_l2 = nullptr) {
    const SpectralVectorField omega_hat = vorticity_hat(grid, velocity_hat);

    RealVectorField velocity = make_real_vector(grid.size());
    RealVectorField omega = make_real_vector(grid.size());
    for (int component = 0; component < 3; ++component) {
        velocity[component] = fft.inverse(velocity_hat[component]);
        omega[component] = fft.inverse(omega_hat[component]);
    }

    RealVectorField rotational = make_real_vector(grid.size());
    for (std::size_t idx = 0; idx < grid.size(); ++idx) {
        rotational[0][idx] =
            velocity[1][idx] * omega[2][idx] - velocity[2][idx] * omega[1][idx];
        rotational[1][idx] =
            velocity[2][idx] * omega[0][idx] - velocity[0][idx] * omega[2][idx];
        rotational[2][idx] =
            velocity[0][idx] * omega[1][idx] - velocity[1][idx] * omega[0][idx];
    }

    SpectralVectorField result = make_spectral_vector(grid.size());
    for (int component = 0; component < 3; ++component) {
        result[component] = fft.forward(rotational[component]);
    }
    project_and_dealias(grid, result);

    if (projected_nonlinear_l2 != nullptr) {
        *projected_nonlinear_l2 = spectral_l2(result);
    }

    for (int i = 0; i < grid.n; ++i) {
        const int kx = grid.wave_number(i);
        for (int j = 0; j < grid.n; ++j) {
            const int ky = grid.wave_number(j);
            for (int k = 0; k < grid.n; ++k) {
                const int kz = grid.wave_number(k);
                const std::size_t idx = grid.index(i, j, k);
                const double k_squared =
                    static_cast<double>(kx * kx + ky * ky + kz * kz);
                for (int component = 0; component < 3; ++component) {
                    result[component][idx] -=
                        viscosity * k_squared * velocity_hat[component][idx];
                }
            }
        }
    }
    return result;
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
    Fft3d& fft,
    const SpectralVectorField& velocity_hat,
    double dt,
    double viscosity) {
    const SpectralVectorField k1 = rhs(grid, fft, velocity_hat, viscosity);
    const SpectralVectorField k2 = rhs(
        grid, fft, add_scaled(velocity_hat, k1, 0.5 * dt), viscosity);
    const SpectralVectorField k3 = rhs(
        grid, fft, add_scaled(velocity_hat, k2, 0.5 * dt), viscosity);
    const SpectralVectorField k4 = rhs(
        grid, fft, add_scaled(velocity_hat, k3, dt), viscosity);

    SpectralVectorField out = velocity_hat;
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

double divergence_l2(
    const Grid& grid,
    const SpectralVectorField& velocity_hat) {
    double sum = 0.0;
    for (int i = 0; i < grid.n; ++i) {
        const int kx = grid.wave_number(i);
        for (int j = 0; j < grid.n; ++j) {
            const int ky = grid.wave_number(j);
            for (int k = 0; k < grid.n; ++k) {
                const int kz = grid.wave_number(k);
                const std::size_t idx = grid.index(i, j, k);
                const Complex divergence =
                    static_cast<double>(kx) * velocity_hat[0][idx]
                    + static_cast<double>(ky) * velocity_hat[1][idx]
                    + static_cast<double>(kz) * velocity_hat[2][idx];
                sum += std::norm(divergence);
            }
        }
    }
    return std::sqrt(sum);
}

double timestep_residual_l2(
    const SpectralVectorField& previous,
    const SpectralVectorField& current,
    const SpectralVectorField& current_rhs,
    double dt) {
    double sum = 0.0;
    for (int component = 0; component < 3; ++component) {
        for (std::size_t i = 0; i < current[component].size(); ++i) {
            const Complex residual =
                (current[component][i] - previous[component][i]) / dt
                - current_rhs[component][i];
            sum += std::norm(residual);
        }
    }
    return std::sqrt(sum);
}

struct Diagnostics {
    double max_velocity{0.0};
    double max_vorticity{0.0};
    double kinetic_energy{0.0};
    double enstrophy{0.0};
    double viscous_dissipation{0.0};
    double divergence_error{0.0};
    double pde_residual{std::numeric_limits<double>::quiet_NaN()};
    double projected_nonlinear{0.0};
    double cfl{0.0};
};

Diagnostics diagnostics(
    const Grid& grid,
    Fft3d& fft,
    const SpectralVectorField& velocity_hat,
    double viscosity,
    double dt,
    double pde_residual) {
    Diagnostics result;
    result.pde_residual = pde_residual;
    result.divergence_error = divergence_l2(grid, velocity_hat);

    const SpectralVectorField omega_hat = vorticity_hat(grid, velocity_hat);
    for (int component = 0; component < 3; ++component) {
        for (const Complex& value : velocity_hat[component]) {
            result.kinetic_energy += 0.5 * std::norm(value);
        }
        for (const Complex& value : omega_hat[component]) {
            result.enstrophy += 0.5 * std::norm(value);
        }
    }
    result.viscous_dissipation = 2.0 * viscosity * result.enstrophy;

    RealVectorField velocity = make_real_vector(grid.size());
    RealVectorField omega = make_real_vector(grid.size());
    for (int component = 0; component < 3; ++component) {
        velocity[component] = fft.inverse(velocity_hat[component]);
        omega[component] = fft.inverse(omega_hat[component]);
    }

    for (std::size_t idx = 0; idx < grid.size(); ++idx) {
        const double speed = std::sqrt(
            velocity[0][idx] * velocity[0][idx]
            + velocity[1][idx] * velocity[1][idx]
            + velocity[2][idx] * velocity[2][idx]);
        const double vorticity_magnitude = std::sqrt(
            omega[0][idx] * omega[0][idx]
            + omega[1][idx] * omega[1][idx]
            + omega[2][idx] * omega[2][idx]);
        result.max_velocity = std::max(result.max_velocity, speed);
        result.max_vorticity = std::max(result.max_vorticity, vorticity_magnitude);
    }

    double projected_nonlinear_l2 = 0.0;
    rhs(grid, fft, velocity_hat, viscosity, &projected_nonlinear_l2);
    result.projected_nonlinear = projected_nonlinear_l2;

    const double dx = 2.0 * std::numbers::pi / static_cast<double>(grid.n);
    result.cfl = result.max_velocity * dt / dx;
    return result;
}

RealVectorField initial_taylor_green(const Grid& grid) {
    RealVectorField initial = make_real_vector(grid.size());
    const double two_pi = 2.0 * std::numbers::pi;

    for (int i = 0; i < grid.n; ++i) {
        const double x = two_pi * static_cast<double>(i) / grid.n;
        for (int j = 0; j < grid.n; ++j) {
            const double y = two_pi * static_cast<double>(j) / grid.n;
            for (int k = 0; k < grid.n; ++k) {
                const std::size_t idx = grid.index(i, j, k);
                initial[0][idx] = std::sin(x) * std::cos(y);
                initial[1][idx] = -std::cos(x) * std::sin(y);
            }
        }
    }
    return initial;
}

void write_diagnostics(
    std::ofstream& out,
    int step,
    double time,
    double dt,
    const Diagnostics& d) {
    out << step << ','
        << time << ','
        << d.max_velocity << ','
        << d.max_vorticity << ','
        << d.kinetic_energy << ','
        << d.enstrophy << ','
        << d.viscous_dissipation << ','
        << d.divergence_error << ','
        << d.pde_residual << ','
        << d.projected_nonlinear << ','
        << dt << ','
        << d.cfl << '\n';
}
}  // namespace

int main(int argc, char** argv) {
    int n = 16;
    double final_time = 0.1;
    double dt = 0.005;
    double viscosity = 0.1;
    std::string output_path = "taylor_green_resolution.csv";

    if (argc > 1) n = std::stoi(argv[1]);
    if (argc > 2) final_time = std::stod(argv[2]);
    if (argc > 3) dt = std::stod(argv[3]);
    if (argc > 4) viscosity = std::stod(argv[4]);
    if (argc > 5) output_path = argv[5];

    if (n < 6 || n % 2 != 0) {
        throw std::invalid_argument("N must be an even integer >= 6");
    }
    if (final_time <= 0.0 || dt <= 0.0 || viscosity <= 0.0) {
        throw std::invalid_argument("final_time, dt, and viscosity must be positive");
    }

    const int steps = static_cast<int>(std::ceil(final_time / dt));
    dt = final_time / static_cast<double>(steps);

    const Grid grid{n};
    Fft3d fft(grid);
    const RealVectorField initial = initial_taylor_green(grid);

    SpectralVectorField velocity_hat = make_spectral_vector(grid.size());
    for (int component = 0; component < 3; ++component) {
        velocity_hat[component] = fft.forward(initial[component]);
    }
    project_and_dealias(grid, velocity_hat);

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("failed to open output file: " + output_path);
    }
    out << std::setprecision(17);
    out << "step,time,max_velocity,max_vorticity,kinetic_energy,enstrophy,"
           "viscous_dissipation,divergence_l2,pde_residual_l2,"
           "projected_nonlinear_l2,dt,cfl\n";

    const Diagnostics initial_diagnostics = diagnostics(
        grid,
        fft,
        velocity_hat,
        viscosity,
        dt,
        std::numeric_limits<double>::quiet_NaN());
    write_diagnostics(out, 0, 0.0, dt, initial_diagnostics);

    for (int step = 1; step <= steps; ++step) {
        const SpectralVectorField previous = velocity_hat;
        velocity_hat = rk4_step(grid, fft, velocity_hat, dt, viscosity);

        const SpectralVectorField current_rhs = rhs(grid, fft, velocity_hat, viscosity);
        const double residual = timestep_residual_l2(
            previous, velocity_hat, current_rhs, dt);
        const Diagnostics d = diagnostics(
            grid, fft, velocity_hat, viscosity, dt, residual);
        write_diagnostics(out, step, step * dt, dt, d);
    }

    RealVectorField numerical = make_real_vector(grid.size());
    for (int component = 0; component < 3; ++component) {
        numerical[component] = fft.inverse(velocity_hat[component]);
    }

    const double exact_factor = std::exp(-2.0 * viscosity * final_time);
    double max_exact_error = 0.0;
    for (int component = 0; component < 3; ++component) {
        for (std::size_t idx = 0; idx < grid.size(); ++idx) {
            max_exact_error = std::max(
                max_exact_error,
                std::abs(numerical[component][idx]
                         - exact_factor * initial[component][idx]));
        }
    }

    const Diagnostics final_diagnostics = diagnostics(
        grid, fft, velocity_hat, viscosity, dt, 0.0);

    std::cout << std::setprecision(16)
              << "FFTW Taylor-Green resolution run\n"
              << "N=" << n << '\n'
              << "steps=" << steps << '\n'
              << "dt=" << dt << '\n'
              << "final_time=" << final_time << '\n'
              << "max_exact_velocity_error=" << max_exact_error << '\n'
              << "final_divergence_l2=" << final_diagnostics.divergence_error << '\n'
              << "final_cfl=" << final_diagnostics.cfl << '\n'
              << "wrote " << output_path << '\n';

    constexpr double regression_tolerance = 1e-8;
    if (!std::isfinite(max_exact_error) ||
        !std::isfinite(final_diagnostics.divergence_error) ||
        max_exact_error > regression_tolerance ||
        final_diagnostics.divergence_error > regression_tolerance) {
        std::cerr << "FFTW Taylor-Green regression failed\n";
        return 1;
    }

    return 0;
}
