// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_traits.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#define NAMESPACE_BEGIN namespace pgm_tpf_hackathon {
#define NAMESPACE_END }

NAMESPACE_BEGIN

using Float = double;
using Dict = std::unordered_map<std::string, Float>;
using Int = int64_t;
using Complex = std::complex<Float>;
using VectorInt = std::vector<Int>;
using VectorScalar = std::vector<Float>;
using VectorComplex = std::vector<Complex>;
using TensorComplex = std::vector<VectorComplex>;

using DenseMatInt = Eigen::Matrix<Int, Eigen::Dynamic, Eigen::Dynamic>;
using DenseMatScalar = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
using DenseMatComplex = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

using SparseMatInt = Eigen::SparseMatrix<Int>;
using SparseMatScalar = Eigen::SparseMatrix<Float>;
using SparseMatComplex = Eigen::SparseMatrix<Complex>;

constexpr Float NaN = std::numeric_limits<Float>::quiet_NaN();
constexpr Float Inf = std::numeric_limits<Float>::infinity();
constexpr Float pi = std::numbers::pi_v<Float>;
constexpr Float BASE_POWER = 1e6;

constexpr Int CONST_POWER = 0;
constexpr Int CONST_CURRENT = 1;
constexpr Int CONST_IMPEDANCE = 2;

using IDx = int64_t; // meh
struct IDx2 {
    IDx x;
    IDx y;
};
struct IDx3 {
    IDx x;
    IDx y;
    IDx z;
};
struct IDx4 {
    IDx w; // meh
    IDx x;
    IDx y;
    IDx z;
};
using VectorIDx = std::vector<IDx>; // meh

struct PgmArray {
    Eigen::MatrixXd data;
    Dict columns;
};

std::ostream& operator<<(std::ostream& os, const PgmArray& pgmArray) {
    os << "Data:\n" << pgmArray.data << "\nColumns:\n";
    for (auto const& [key, value] : pgmArray.columns) {
        os << key << ": " << value << "\n";
    }
    return os;
}

using PgmData = std::unordered_map<std::string, PgmArray>;
using PgmDataset = std::unordered_map<std::string, PgmData>;
using PgmBatchDataset = std::unordered_map<std::string, std::vector<PgmDataset>>;
using PgmResultType = std::unordered_map<std::string, PgmData>;

NAMESPACE_END
