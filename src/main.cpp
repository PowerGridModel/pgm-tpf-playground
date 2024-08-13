// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <common.hpp>
#include <utils.hpp>

#include <tpf.hpp>

int main() {
    auto const dataset = pgm_tpf_hackathon::generate_fictional_grid(3, 5, 0.1, 1.0, 1000, 500, 0.95, 10, 0.8, 1.2);
    auto const pgm_data = dataset.at("pgm_data");
    auto const update_data = dataset.at("pgm_update_data");

    // Accessing elements from pgm_data
    try {
        auto some_value = pgm_data.at("line");
        std::cout << "Value: \n" << some_value << std::endl;
    } catch (std::out_of_range const& e) {
        std::cerr << "Key not found: " << e.what() << std::endl;
    }

    // Accessing elements from update_data
    try {
        auto another_value_id = update_data.at("id");
        auto another_value_ps = update_data.at("p_specified");
        auto another_value_qs = update_data.at("q_specified");
        std::cout << "Value: \n" << another_value_id << std::endl;
        std::cout << "Value: \n" << another_value_ps << std::endl;
        std::cout << "Value: \n" << another_value_qs << std::endl;
    } catch (std::out_of_range const& e) {
        std::cerr << "Key not found: " << e.what() << std::endl;
    }

    pgm_tpf_hackathon::TPF tpf{pgm_data, 50.0};

    auto const result = tpf.calculate_power_flow(update_data);

    return 0;
}
