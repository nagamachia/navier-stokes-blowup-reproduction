#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    double h = 0.005;
    int decades = 8;
    int samples_per_decade = 20;
    std::string output_path = "core_scaling.csv";

    if (argc > 1) h = std::stod(argv[1]);
    if (argc > 2) decades = std::stoi(argv[2]);
    if (argc > 3) samples_per_decade = std::stoi(argv[3]);
    if (argc > 4) output_path = argv[4];

    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }
    if (decades <= 0 || samples_per_decade <= 0) {
        throw std::invalid_argument("decades and samples_per_decade must be positive");
    }

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("failed to open output file: " + output_path);
    }

    out << std::setprecision(17);
    out << "tau,t,L_r,L_z,U,U_r_scale,V_core,E_core,E_reference,relative_error\n";

    const int total_samples = decades * samples_per_decade + 1;
    double max_relative_error = 0.0;

    for (int i = 0; i < total_samples; ++i) {
        const double log10_tau = -static_cast<double>(i) / samples_per_decade;
        const double tau = std::pow(10.0, log10_tau);
        const double t = 1.0 - tau;

        const double L_r = std::pow(tau, 0.5);
        const double L_z = std::pow(tau, 0.5 - h);
        const double U = std::pow(tau, -0.5 - h);
        const double U_r_scale = std::pow(tau, -0.5);

        const double V_core = L_r * L_r * L_z;
        const double E_core = U * U * V_core;
        const double E_reference = std::pow(tau, 0.5 - 3.0 * h);
        const double relative_error = std::abs(E_core - E_reference) / E_reference;

        if (relative_error > max_relative_error) {
            max_relative_error = relative_error;
        }

        out << tau << ','
            << t << ','
            << L_r << ','
            << L_z << ','
            << U << ','
            << U_r_scale << ','
            << V_core << ','
            << E_core << ','
            << E_reference << ','
            << relative_error << '\n';
    }

    std::cout << std::setprecision(10);
    std::cout << "Phase 0 core-scaling experiment\n";
    std::cout << "h = " << h << '\n';
    std::cout << "velocity exponent = " << (-0.5 - h) << '\n';
    std::cout << "energy exponent = " << (0.5 - 3.0 * h) << '\n';
    std::cout << "samples = " << total_samples << '\n';
    std::cout << "max energy identity relative error = " << max_relative_error << '\n';
    std::cout << "wrote " << output_path << '\n';

    return 0;
}
