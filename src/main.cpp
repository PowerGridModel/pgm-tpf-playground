// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <common.hpp>
#include <utils.hpp>

#include <tpf.hpp>

#include <ctime>

// Example: PGM_TPF_Hackathon_2024.[exe/o] 3 5 0.1 1.0 1000 500 0.95 10 0.8 1.2 0 0
int main(int argc, char* argv[]) {
    auto const& options = pgm_tpf_hackathon::parse_experiment_options(argc, argv);
    auto const& dataset = pgm_tpf_hackathon::generate_fictional_grid(options);
    auto const& pgm_data = dataset.at("pgm_data");
    auto const& update_data = dataset.at("pgm_update_data");

    print_pgm_dataset(dataset, options.print_res);

    clock_t start_init = clock();
    pgm_tpf_hackathon::TPF tpf{pgm_data, pgm_tpf_hackathon::frequency};
    clock_t end_init = clock();
    double init_duration = 1000000.0 * (end_init - start_init) / CLOCKS_PER_SEC;
    std::cout << "TPF initialization took " << init_duration << " us\n";

    clock_t start_calc = clock();
    auto const result = tpf.calculate_power_flow(update_data);
    clock_t end_calc = clock();
    double calc_duration = 1000000.0 * (end_calc - start_calc) / CLOCKS_PER_SEC;
    std::cout << "TPF power flow calculation took " << calc_duration << " us\n";

    print_tpf_result(result, options.print_res);

    return 0;
}
