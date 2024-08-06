#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#define NAMESPACE_BEGIN namespace pgm_tpf_hackathon {
#define NAMESPACE_END }

NAMESPACE_BEGIN

using Float = double;
using Dict = std::unordered_map<std::string, Float>;
using Int = int64_t;
using VectorInt = std::vector<Int>;
using VectorScalar = std::vector<Float>;
using VectorComplex = std::vector<std::complex<Float>>;
using TensorComplex = std::vector<VectorComplex>;

#define NaN std::nan("")
#define Inf std::numeric_limits<Float>::infinity()
constexpr Float M_PI = std::numbers::pi_v<Float>;
#define CONST_POWER 0
#define CONST_CURRENT 1
#define CONST_IMPEDANCE 2
#define BASE_POWER 1e6

using IDx = int64_t;
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
    IDx w;  // meh
    IDx x;
    IDx y;
    IDx z;
};

struct PgmArray {
    Eigen::MatrixXd data;
    Dict columns;
};

std::ostream &operator<<(std::ostream &os, const PgmArray &pgmArray) {
    os << "Data:\n" << pgmArray.data << "\nColumns:\n";
    for (const auto &[key, value] : pgmArray.columns) {
        os << key << ": " << value << "\n";
    }
    return os;
}

using PgmData = std::unordered_map<std::string, PgmArray>;
using PgmDataset = std::unordered_map<std::string, PgmData>;

NAMESPACE_END
