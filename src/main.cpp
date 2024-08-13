// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <common.hpp>
#include <utils.hpp>

#include <tpf.hpp>

#include <chrono>

// Example: PGM_TPF_Hackathon_2024.[exe] 3 5 0.1 1.0 1000 500 0.95 10 0.8 1.2
int main(int argc, char* argv[]) {
    auto const& options = pgm_tpf_hackathon::parse_experiment_options(argc, argv);
    auto const& dataset = pgm_tpf_hackathon::generate_fictional_grid(options);
    auto const& pgm_data = dataset.at("pgm_data");
    auto const& update_data = dataset.at("pgm_update_data");

    print_pgm_dataset(dataset);

    auto start_init = std::chrono::high_resolution_clock::now();
    pgm_tpf_hackathon::TPF tpf{pgm_data, pgm_tpf_hackathon::frequency};
    auto end_init = std::chrono::high_resolution_clock::now();
    auto init_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_init - start_init).count();
    std::cout << "TPF initialization took " << init_duration << " us\n";

    auto start_calc = std::chrono::high_resolution_clock::now();
    auto const result = tpf.calculate_power_flow(update_data);
    auto end_calc = std::chrono::high_resolution_clock::now();
    auto calc_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_calc - start_calc).count();
    std::cout << "TPF power flow calculation took " << calc_duration << " us\n";

    print_tpf_result(result);

    return 0;
}
