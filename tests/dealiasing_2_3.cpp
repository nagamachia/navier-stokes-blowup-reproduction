#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

namespace {
using Complex = std::complex<double>;

int wave_number(int index, int n) {
    return (index <= n / 2) ? index : index - n;
}

std::vector<Complex> dft(const std::vector<double>& values) {
    const int n = static_cast<int>(values.size());
    const double two_pi = 2.0 * std::numbers::pi;
    std::vector<Complex> spectrum(n, {0.0, 0.0});
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            const double angle = -two_pi * static_cast<double>(k * j) / static_cast<double>(n);
            spectrum[k] += values[j] * Complex{std::cos(angle), std::sin(angle)};
        }
        spectrum[k] /= static_cast<double>(n);
    }
    return spectrum;
}

void apply_two_thirds_filter(std::vector<Complex>& spectrum) {
    const int n = static_cast<int>(spectrum.size());
    const int cutoff = n / 3;
    for (int k = 0; k < n; ++k) {
        if (std::abs(wave_number(k, n)) > cutoff) {
            spectrum[k] = {0.0, 0.0};
        }
    }
}
}  // namespace

int main() {
    constexpr int n = 12;
    constexpr double tolerance = 1e-12;
    const double two_pi = 2.0 * std::numbers::pi;

    std::vector<double> product(n, 0.0);
    for (int j = 0; j < n; ++j) {
        const double x = two_pi * static_cast<double>(j) / static_cast<double>(n);
        product[j] = 2.0 * std::cos(3.0 * x) * std::cos(4.0 * x);
    }

    auto spectrum = dft(product);
    const int alias_index = n - 5;  // k = -5, alias of the unresolved k = +7 mode.
    const double alias_before = std::abs(spectrum[alias_index]);

    apply_two_thirds_filter(spectrum);

    double max_unexpected = 0.0;
    for (int k = 0; k < n; ++k) {
        const int wn = wave_number(k, n);
        if (std::abs(wn) == 1) {
            max_unexpected = std::max(max_unexpected, std::abs(spectrum[k] - Complex{0.5, 0.0}));
        } else {
            max_unexpected = std::max(max_unexpected, std::abs(spectrum[k]));
        }
    }

    std::cout << std::setprecision(16)
              << "N=" << n << '\n'
              << "cutoff=" << n / 3 << '\n'
              << "alias_amplitude_before_filter=" << alias_before << '\n'
              << "alias_amplitude_after_filter=" << std::abs(spectrum[alias_index]) << '\n'
              << "max_filtered_spectrum_error=" << max_unexpected << '\n';

    if (std::abs(alias_before - 0.5) > tolerance || max_unexpected > tolerance) {
        std::cerr << "2/3 dealiasing check failed\n";
        return 1;
    }

    std::cout << "2/3 dealiasing check passed\n";
    return 0;
}
