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
        auto some_value = pgm_data.at("source");
        std::cout << "Value: \n" << some_value << std::endl;
    } catch (std::out_of_range const& e) {
        std::cerr << "Key not found: " << e.what() << std::endl;
    }

    // Accessing elements from update_data
    try {
        auto another_value = update_data.at("sym_load");
        std::cout << "Value: \n" << another_value << std::endl;
    } catch (std::out_of_range const& e) {
        std::cerr << "Key not found: " << e.what() << std::endl;
    }

    return 0;
}
