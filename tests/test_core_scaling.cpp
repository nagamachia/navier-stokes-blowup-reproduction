#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool approx(double a, double b, double rel_tol = 1e-12) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel_tol * scale;
}

}

int main() {
    const double h = 0.005;
    const double tau = 1e-4;

    const double lr = std::pow(tau, 0.5);
    const double lz = std::pow(tau, 0.5 - h);
    const double u = std::pow(tau, -0.5 - h);
    const double volume = lr * lr * lz;
    const double energy = u * u * volume;

    assert(approx(lr, std::pow(tau, 0.5)));
    assert(approx(lz, std::pow(tau, 0.5 - h)));
    assert(approx(u, std::pow(tau, -0.5 - h)));
    assert(approx(energy, std::pow(tau, 0.5 - 3.0 * h)));

    const double tau2 = 1e-6;
    const double u2 = std::pow(tau2, -0.5 - h);
    const double e2 = std::pow(tau2, 0.5 - 3.0 * h);

    assert(u2 > u);
    assert(e2 < energy);

    std::cout << "core scaling regression checks passed\n";
    return 0;
}
