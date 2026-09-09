#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool approx(double a, double b, double rel_tol = 1e-12) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= rel_tol * scale;
}

bool require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "core scaling regression failed: " << message << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main() {
    const double h = 0.005;
    const double tau = 1e-4;

    const double lr = std::pow(tau, 0.5);
    const double lz = std::pow(tau, 0.5 - h);
    const double u = std::pow(tau, -0.5 - h);
    const double volume = lr * lr * lz;
    const double energy = u * u * volume;

    bool ok = true;
    ok = require(approx(lr, std::pow(tau, 0.5)), "radial length scale") && ok;
    ok = require(approx(lz, std::pow(tau, 0.5 - h)), "axial length scale") && ok;
    ok = require(approx(u, std::pow(tau, -0.5 - h)), "velocity scale") && ok;
    ok = require(
        approx(energy, std::pow(tau, 0.5 - 3.0 * h)),
        "core energy power law") && ok;

    const double tau2 = 1e-6;
    const double u2 = std::pow(tau2, -0.5 - h);
    const double e2 = std::pow(tau2, 0.5 - 3.0 * h);

    ok = require(u2 > u, "velocity should grow as tau decreases") && ok;
    ok = require(e2 < energy, "core energy should decrease as tau decreases") && ok;

    if (!ok) {
        return 1;
    }

    std::cout << "core scaling regression checks passed\n";
    return 0;
}
