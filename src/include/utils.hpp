// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "./common.hpp"

NAMESPACE_BEGIN

// Fictional Grid Generator

// Indices
int const COL_R1 = 5;
int const COL_X1 = 6;
int const COL_C1 = 7;
int const COL_C0 = 8;
int const COL_TAN1 = 9;
int const COL_TAN0 = 10;
int const COL_I_N = 11;

// Constants
Float const u_rated = 10e3;
Float const frequency = 50.0;

// Source parameters
Float const source_sk = 1e20;
Float const source_rx = 0.1;
Float const source_01 = 1.0;
Float const source_u_ref = 1.0;
int const source_node = 0;

// Cable parameters per km
std::string const cable_type = "630Al";
Dict const cable_param = {{"r1", 0.063},   {"x1", 0.103},   {"c1", 0.5e-6}, {"c0", 0.3e-6},
                          {"tan1", 0.003}, {"tan0", 0.001}, {"i_n", NaN}};

Dict const cable_param_pp = {{"c_nf_per_km", cable_param.at("c1") * 1e9},
                             {"r_ohm_per_km", cable_param.at("r1")},
                             {"x_ohm_per_km", cable_param.at("x1")},
                             {"g_us_per_km", cable_param.at("tan1") * cable_param.at("c1") * 2 * pi* frequency * 1e6},
                             {"c0_nf_per_km", cable_param.at("c0") * 1e9},
                             {"g0_us_per_km", cable_param.at("tan0") * cable_param.at("c0") * 2 * pi* frequency * 1e6},
                             {"max_i_ka", cable_param.at("i_n") * 1e-3}};

PgmArray initialize_array(int rows, int cols = 0) {
    PgmArray array;
    array.data = Eigen::MatrixXd::Zero(rows, cols);
    return array;
}

PgmDataset generate_fictional_grid(int n_feeder, int n_node_per_feeder, Float cable_length_km_min,
                                   Float cable_length_km_max, Float load_p_w_max, Float load_p_w_min, Float pf,
                                   int n_step, Float load_scaling_min, Float load_scaling_max, int seed = 0) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Float> dist_length(cable_length_km_min, cable_length_km_max);
    std::uniform_real_distribution<Float> dist_load(load_p_w_min / 3.0, load_p_w_max / 3.0);
    std::uniform_real_distribution<Float> dist_scaling(load_scaling_min, load_scaling_max);

    int n_node = n_feeder * n_node_per_feeder + 1;
    PgmData pgm_data;

    // Node
    pgm_data["node"] = initialize_array(n_node, 2);
    pgm_data["node"].data.col(0) = Eigen::VectorXd::LinSpaced(n_node, 0, n_node - 1);
    pgm_data["node"].data.col(1).setConstant(u_rated);

    // Line
    int n_line = n_node - 1;

    auto to_node_raw = [&n_node_per_feeder, &n_feeder]() {
        VectorInt to_node_feeder(n_node_per_feeder);
        std::iota(to_node_feeder.begin(), to_node_feeder.end(), 1);

        std::vector<VectorInt> reshaped_to_node_feeder(n_feeder, VectorInt(n_node_per_feeder));
        for (int i = 0; i < n_feeder; ++i) {
            for (int j = 0; j < n_node_per_feeder; ++j) {
                reshaped_to_node_feeder[i][j] = to_node_feeder[j] + i * n_node_per_feeder;
            }
        }

        VectorReal to_nodes;
        for (const auto& feeder : reshaped_to_node_feeder) {
            to_nodes.insert(to_nodes.end(), feeder.begin(), feeder.end());
        }

        return to_nodes;
    }();
    Eigen::VectorXd to_node = Eigen::Map<Eigen::VectorXd>(to_node_raw.data(), to_node_raw.size());

    auto from_node_raw = [&n_node_per_feeder, &n_feeder]() {
        VectorInt from_node_feeder(n_node_per_feeder - 1);
        std::iota(from_node_feeder.begin(), from_node_feeder.end(), 1);

        std::vector<VectorInt> reshaped_from_node_feeder(n_feeder, VectorInt(n_node_per_feeder - 1));
        for (int i = 0; i < n_feeder; ++i) {
            for (int j = 0; j < n_node_per_feeder - 1; ++j) {
                reshaped_from_node_feeder[i][j] = from_node_feeder[j] + i * n_node_per_feeder;
            }
        }

        for (int i = 0; i < n_feeder; ++i) {
            reshaped_from_node_feeder[i].insert(reshaped_from_node_feeder[i].begin(), 0);
        }

        VectorReal from_nodes;
        for (auto const& feeder : reshaped_from_node_feeder) {
            from_nodes.insert(from_nodes.end(), feeder.begin(), feeder.end());
        }

        return from_nodes;
    }();
    Eigen::VectorXd from_node = Eigen::Map<Eigen::VectorXd>(from_node_raw.data(), from_node_raw.size());

    Eigen::VectorXd length = Eigen::VectorXd::NullaryExpr(n_line, [&]() { return dist_length(rng); });

    pgm_data["line"] = initialize_array(n_line, 12);
    pgm_data["line"].data.col(0) = Eigen::VectorXd::LinSpaced(n_line, n_node, n_node + n_line - 1);
    pgm_data["line"].data.col(1) = from_node;
    pgm_data["line"].data.col(2) = to_node;
    pgm_data["line"].data.col(3).setConstant(1);
    pgm_data["line"].data.col(4).setConstant(1);

    pgm_data["line"].data.col(COL_R1) = Eigen::VectorXd::Constant(n_line, cable_param.at("r1")).cwiseProduct(length);
    pgm_data["line"].data.col(COL_X1) = Eigen::VectorXd::Constant(n_line, cable_param.at("x1")).cwiseProduct(length);
    pgm_data["line"].data.col(COL_C1) = Eigen::VectorXd::Constant(n_line, cable_param.at("c1")).cwiseProduct(length);
    pgm_data["line"].data.col(COL_C0) = Eigen::VectorXd::Constant(n_line, cable_param.at("c0")).cwiseProduct(length);
    pgm_data["line"].data.col(COL_TAN1).setConstant(cable_param.at("tan1"));
    pgm_data["line"].data.col(COL_TAN0).setConstant(cable_param.at("tan0"));
    pgm_data["line"].data.col(COL_I_N).setConstant(cable_param.at("i_n"));

    // Load
    int n_load = n_node - 1;
    pgm_data["sym_load"] = initialize_array(n_load, 6);
    pgm_data["sym_load"].data.col(0) =
        Eigen::VectorXd::LinSpaced(n_load, n_node + n_line, n_node + n_line + n_load - 1);
    pgm_data["sym_load"].data.col(1) = pgm_data["node"].data.col(0).tail(n_load);
    pgm_data["sym_load"].data.col(2).setConstant(1); // status
    pgm_data["sym_load"].data.col(3).setConstant(1); // Assuming const_power type
    pgm_data["sym_load"].data.col(4) = Eigen::VectorXd::NullaryExpr(n_load, [&]() { return dist_load(rng); });
    pgm_data["sym_load"].data.col(5) =
        pgm_data["sym_load"].data.col(4).cwiseProduct(Eigen::VectorXd::Constant(n_load, std::sqrt(1 - pf * pf) / pf));

    // Source
    int source_id = n_node + n_line + n_load;
    pgm_data["source"] = initialize_array(1, 7);
    pgm_data["source"].data(0, 0) = source_id;
    pgm_data["source"].data(0, 1) = source_node;
    pgm_data["source"].data(0, 2) = 1;
    pgm_data["source"].data(0, 3) = source_u_ref;
    pgm_data["source"].data(0, 4) = source_sk;
    pgm_data["source"].data(0, 5) = source_rx;
    pgm_data["source"].data(0, 6) = source_01;

    // Generate time series
    Eigen::MatrixXd scaling = Eigen::MatrixXd::NullaryExpr(n_step, n_load, [&]() { return dist_scaling(rng); });

    PgmArray sym_load_profile_id = initialize_array(n_step, n_load);
    PgmArray sym_load_profile_p_specified = initialize_array(n_step, n_load);
    PgmArray sym_load_profile_q_specified = initialize_array(n_step, n_load);

    sym_load_profile_id.data = pgm_data["sym_load"].data.col(0).transpose().replicate(n_step, 1);
    sym_load_profile_p_specified.data =
        pgm_data["sym_load"].data.col(4).transpose().replicate(n_step, 1).cwiseProduct(scaling);
    sym_load_profile_q_specified.data =
        pgm_data["sym_load"].data.col(5).transpose().replicate(n_step, 1).cwiseProduct(scaling);

    // "update" "sym_load" are omitted
    PgmData pgm_update_dataset = {{{"id", sym_load_profile_id},
                                   {"p_specified", sym_load_profile_p_specified},
                                   {"q_specified", sym_load_profile_q_specified}}};

    return PgmDataset{{"pgm_data", pgm_data}, {"pgm_update_data", pgm_update_dataset}};
}

NAMESPACE_END
