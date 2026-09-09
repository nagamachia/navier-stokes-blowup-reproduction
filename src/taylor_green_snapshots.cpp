#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;
using Field = std::vector<Complex>;

struct Grid {
    int n;
    std::size_t size() const { return static_cast<std::size_t>(n) * n; }
    std::size_t idx(int i, int j) const { return static_cast<std::size_t>(i) * n + j; }
    int k(int i) const { return (i <= n / 2) ? i : i - n; }
};

class Fft2d {
public:
    explicit Fft2d(const Grid& g) : g_(g) {
        in_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * g_.size()));
        out_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * g_.size()));
        if (!in_ || !out_) throw std::runtime_error("FFTW allocation failed");
        fwd_ = fftw_plan_dft_2d(g_.n, g_.n, in_, out_, FFTW_FORWARD, FFTW_ESTIMATE);
        inv_ = fftw_plan_dft_2d(g_.n, g_.n, in_, out_, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (!fwd_ || !inv_) throw std::runtime_error("FFTW plan creation failed");
    }
    ~Fft2d() {
        if (fwd_) fftw_destroy_plan(fwd_);
        if (inv_) fftw_destroy_plan(inv_);
        fftw_free(in_);
        fftw_free(out_);
    }

    Field forward_real(const std::vector<double>& a) {
        for (std::size_t p = 0; p < g_.size(); ++p) {
            in_[p][0] = a[p]; in_[p][1] = 0.0;
        }
        fftw_execute(fwd_);
        Field h(g_.size());
        const double s = 1.0 / static_cast<double>(g_.size());
        for (std::size_t p = 0; p < g_.size(); ++p) h[p] = s * Complex(out_[p][0], out_[p][1]);
        return h;
    }

    std::vector<double> inverse_real(const Field& h) {
        for (std::size_t p = 0; p < g_.size(); ++p) {
            in_[p][0] = h[p].real(); in_[p][1] = h[p].imag();
        }
        fftw_execute(inv_);
        std::vector<double> a(g_.size());
        for (std::size_t p = 0; p < g_.size(); ++p) a[p] = out_[p][0];
        return a;
    }

private:
    Grid g_;
    fftw_complex* in_{nullptr};
    fftw_complex* out_{nullptr};
    fftw_plan fwd_{nullptr};
    fftw_plan inv_{nullptr};
};

struct State { Field u, v; };

void project_dealias(const Grid& g, State& s) {
    const int cutoff = g.n / 3;
    for (int i = 0; i < g.n; ++i) {
        const int kx = g.k(i);
        for (int j = 0; j < g.n; ++j) {
            const int ky = g.k(j);
            const auto p = g.idx(i, j);
            if (std::abs(kx) > cutoff || std::abs(ky) > cutoff) {
                s.u[p] = 0.0; s.v[p] = 0.0; continue;
            }
            const double k2 = static_cast<double>(kx * kx + ky * ky);
            if (k2 == 0.0) continue;
            const Complex dot = static_cast<double>(kx) * s.u[p] + static_cast<double>(ky) * s.v[p];
            s.u[p] -= static_cast<double>(kx) * dot / k2;
            s.v[p] -= static_cast<double>(ky) * dot / k2;
        }
    }
}

State rhs(const Grid& g, Fft2d& fft, const State& s, double nu) {
    const Complex I{0.0, 1.0};
    Field omega_hat(g.size());
    for (int i = 0; i < g.n; ++i) {
        const double kx = g.k(i);
        for (int j = 0; j < g.n; ++j) {
            const double ky = g.k(j);
            const auto p = g.idx(i, j);
            omega_hat[p] = I * (kx * s.v[p] - ky * s.u[p]);
        }
    }
    const auto u = fft.inverse_real(s.u);
    const auto v = fft.inverse_real(s.v);
    const auto omega = fft.inverse_real(omega_hat);
    std::vector<double> nx(g.size()), ny(g.size());
    for (std::size_t p = 0; p < g.size(); ++p) {
        nx[p] = v[p] * omega[p];
        ny[p] = -u[p] * omega[p];
    }
    State r{fft.forward_real(nx), fft.forward_real(ny)};
    project_dealias(g, r);
    for (int i = 0; i < g.n; ++i) {
        const int kx = g.k(i);
        for (int j = 0; j < g.n; ++j) {
            const int ky = g.k(j);
            const auto p = g.idx(i, j);
            const double k2 = static_cast<double>(kx * kx + ky * ky);
            r.u[p] -= nu * k2 * s.u[p];
            r.v[p] -= nu * k2 * s.v[p];
        }
    }
    return r;
}

State add_scaled(const State& a, const State& b, double scale) {
    State r = a;
    for (std::size_t p = 0; p < a.u.size(); ++p) {
        r.u[p] += scale * b.u[p];
        r.v[p] += scale * b.v[p];
    }
    return r;
}

State rk4(const Grid& g, Fft2d& fft, const State& s, double dt, double nu) {
    const State k1 = rhs(g, fft, s, nu);
    const State k2 = rhs(g, fft, add_scaled(s, k1, 0.5 * dt), nu);
    const State k3 = rhs(g, fft, add_scaled(s, k2, 0.5 * dt), nu);
    const State k4 = rhs(g, fft, add_scaled(s, k3, dt), nu);
    State r = s;
    for (std::size_t p = 0; p < s.u.size(); ++p) {
        r.u[p] += dt * (k1.u[p] + 2.0 * k2.u[p] + 2.0 * k3.u[p] + k4.u[p]) / 6.0;
        r.v[p] += dt * (k1.v[p] + 2.0 * k2.v[p] + 2.0 * k3.v[p] + k4.v[p]) / 6.0;
    }
    project_dealias(g, r);
    return r;
}

void write_snapshot(std::ofstream& out, const Grid& g, Fft2d& fft, const State& s,
                    double time, double nu) {
    const Complex I{0.0, 1.0};
    Field omega_hat(g.size());
    for (int i = 0; i < g.n; ++i) {
        const double kx = g.k(i);
        for (int j = 0; j < g.n; ++j) {
            const double ky = g.k(j);
            const auto p = g.idx(i, j);
            omega_hat[p] = I * (kx * s.v[p] - ky * s.u[p]);
        }
    }
    const auto u = fft.inverse_real(s.u);
    const auto v = fft.inverse_real(s.v);
    const auto omega = fft.inverse_real(omega_hat);
    const double a = std::exp(-2.0 * nu * time);
    const double two_pi = 2.0 * std::numbers::pi;
    for (int i = 0; i < g.n; ++i) {
        const double x = two_pi * static_cast<double>(i) / g.n;
        for (int j = 0; j < g.n; ++j) {
            const double y = two_pi * static_cast<double>(j) / g.n;
            const auto p = g.idx(i, j);
            const double ue = a * std::sin(x) * std::cos(y);
            const double ve = -a * std::cos(x) * std::sin(y);
            const double err = std::hypot(u[p] - ue, v[p] - ve);
            out << time << ',' << i << ',' << j << ',' << x << ',' << y << ','
                << u[p] << ',' << v[p] << ',' << omega[p] << ','
                << ue << ',' << ve << ',' << err << '\n';
        }
    }
}
}  // namespace

int main(int argc, char** argv) {
    int n = 32;
    double final_time = 0.1;
    double dt = 0.0025;
    double nu = 0.1;
    std::string output = "taylor_green_numerical_snapshots.csv";
    if (argc > 1) n = std::stoi(argv[1]);
    if (argc > 2) final_time = std::stod(argv[2]);
    if (argc > 3) dt = std::stod(argv[3]);
    if (argc > 4) nu = std::stod(argv[4]);
    if (argc > 5) output = argv[5];
    if (n < 8 || n % 2 != 0) throw std::invalid_argument("N must be even and >= 8");

    Grid g{n};
    Fft2d fft(g);
    std::vector<double> u0(g.size()), v0(g.size());
    const double two_pi = 2.0 * std::numbers::pi;
    for (int i = 0; i < n; ++i) {
        const double x = two_pi * static_cast<double>(i) / n;
        for (int j = 0; j < n; ++j) {
            const double y = two_pi * static_cast<double>(j) / n;
            const auto p = g.idx(i, j);
            u0[p] = std::sin(x) * std::cos(y);
            v0[p] = -std::cos(x) * std::sin(y);
        }
    }
    State state{fft.forward_real(u0), fft.forward_real(v0)};
    project_dealias(g, state);

    std::ofstream out(output);
    if (!out) throw std::runtime_error("failed to open output");
    out << std::setprecision(17);
    out << "time,i,j,x,y,u,v,omega_z,exact_u,exact_v,velocity_error\n";

    const std::vector<double> targets{0.0, 0.25 * final_time, 0.5 * final_time,
                                      0.75 * final_time, final_time};
    std::size_t next_target = 0;
    double t = 0.0;
    write_snapshot(out, g, fft, state, t, nu);
    next_target = 1;

    while (t < final_time - 1e-15) {
        const double step_dt = std::min(dt, final_time - t);
        state = rk4(g, fft, state, step_dt, nu);
        t += step_dt;
        while (next_target < targets.size() && t + 1e-12 >= targets[next_target]) {
            write_snapshot(out, g, fft, state, t, nu);
            ++next_target;
        }
    }

    std::cout << "wrote numerical Taylor-Green snapshots: " << output << '\n';
    return 0;
}
