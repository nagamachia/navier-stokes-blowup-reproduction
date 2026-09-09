#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

int main() {
    constexpr int n = 32;
    constexpr double two_pi = 2.0 * std::numbers::pi;

    std::vector<double> x(n);
    std::vector<double> f(n);
    std::vector<std::complex<double>> fhat(n, {0.0, 0.0});
    std::vector<std::complex<double>> dfhat(n, {0.0, 0.0});
    std::vector<double> df(n, 0.0);

    for (int j = 0; j < n; ++j) {
        x[j] = two_pi * static_cast<double>(j) / static_cast<double>(n);
        f[j] = std::sin(3.0 * x[j]) + 0.5 * std::cos(5.0 * x[j]);
    }

    // Reference O(N^2) DFT. This is intentionally simple and auditable.
    // It will later be replaced by FFTW without changing the spectral rule.
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            const double angle = -two_pi * static_cast<double>(k * j) / static_cast<double>(n);
            fhat[k] += f[j] * std::complex<double>(std::cos(angle), std::sin(angle));
        }
        fhat[k] /= static_cast<double>(n);
    }

    for (int k = 0; k < n; ++k) {
        const int wave_number = (k <= n / 2) ? k : k - n;
        dfhat[k] = std::complex<double>(0.0, static_cast<double>(wave_number)) * fhat[k];
    }

    for (int j = 0; j < n; ++j) {
        std::complex<double> value{0.0, 0.0};
        for (int k = 0; k < n; ++k) {
            const double angle = two_pi * static_cast<double>(k * j) / static_cast<double>(n);
            value += dfhat[k] * std::complex<double>(std::cos(angle), std::sin(angle));
        }
        df[j] = value.real();
    }

    double max_error = 0.0;
    for (int j = 0; j < n; ++j) {
        const double exact = 3.0 * std::cos(3.0 * x[j]) - 2.5 * std::sin(5.0 * x[j]);
        max_error = std::max(max_error, std::abs(df[j] - exact));
    }

    std::cout << std::setprecision(16)
              << "N=" << n << "\n"
              << "max_abs_error=" << max_error << "\n";

    constexpr double tolerance = 1e-11;
    if (max_error > tolerance) {
        std::cerr << "Fourier differentiation check failed\n";
        return 1;
    }

    std::cout << "Fourier differentiation check passed\n";
    return 0;
}
