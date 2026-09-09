#include "similarity_coordinates.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    double h = 0.005;
    int decades = 6;
    int eta_samples = 19;
    std::string output_path = "similarity_coordinates.csv";

    if (argc > 1) h = std::stod(argv[1]);
    if (argc > 2) decades = std::stoi(argv[2]);
    if (argc > 3) eta_samples = std::stoi(argv[3]);
    if (argc > 4) output_path = argv[4];

    if (!(h > 0.0 && h < 0.01)) {
        throw std::invalid_argument("h must satisfy 0 < h < 0.01");
    }
    if (decades <= 0 || eta_samples < 3) {
        throw std::invalid_argument("decades must be positive and eta_samples >= 3");
    }

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("failed to open output file: " + output_path);
    }
    out << std::setprecision(17);
    out << "tau,eta_target,z,q_reference,q_solved,eta_solved,relative_q_error,eta_error\n";

    const double D = 0.5 - h;
    double max_q_error = 0.0;
    double max_eta_error = 0.0;

    for (int decade = 0; decade <= decades; ++decade) {
        const double tau = std::pow(10.0, -static_cast<double>(decade));
        for (int j = 0; j < eta_samples; ++j) {
            const double eta = -0.9 + 1.8 * static_cast<double>(j) /
                static_cast<double>(eta_samples - 1);
            const double q_reference = tau / (1.0 - eta * eta);
            const double z = std::pow(q_reference, D) * eta;
            const double q_solved = nsblowup::solve_q(tau, z, h);
            const double eta_solved = z / std::pow(q_solved, D);
            const double q_error = std::abs(q_solved - q_reference) / q_reference;
            const double eta_error = std::abs(eta_solved - eta);
            max_q_error = std::max(max_q_error, q_error);
            max_eta_error = std::max(max_eta_error, eta_error);

            out << tau << ',' << eta << ',' << z << ',' << q_reference << ','
                << q_solved << ',' << eta_solved << ',' << q_error << ','
                << eta_error << '\n';
        }
    }

    std::cout << std::setprecision(12)
              << "similarity coordinate diagnostic\n"
              << "h=" << h << '\n'
              << "max_relative_q_error=" << max_q_error << '\n'
              << "max_eta_error=" << max_eta_error << '\n'
              << "wrote " << output_path << '\n';

    if (max_q_error > 5e-13 || max_eta_error > 5e-13) {
        std::cerr << "similarity coordinate regression failed\n";
        return 1;
    }
    return 0;
}
