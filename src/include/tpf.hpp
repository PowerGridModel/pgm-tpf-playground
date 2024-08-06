#pragma once

#include "./common.hpp"
#include "./utils.hpp"

NAMESPACE_BEGIN

// Minimal TPF impl ==============================

void set_load_pu(std::vector<std::complex<double>> &load_pu,
                 const std::vector<double> &p_array,
                 const std::vector<double> &q_array) {
    std::transform(p_array.begin(), p_array.end(), q_array.begin(),
                   load_pu.begin(), [](double p, double q) {
                       return std::complex<double>(p, q) / BASE_POWER;
                   });
}

std::vector<std::complex<double>> get_load_pu(
    const std::vector<double> &p_specified,
    const std::vector<double> &q_specified) {
    std::vector<std::complex<double>> load_pu(p_specified.size());
    set_load_pu(load_pu, p_specified, q_specified);
    return load_pu;
}

void set_rhs_impl(std::vector<std::vector<std::complex<double>>> &rhs,
                  const std::vector<std::complex<double>> &load_pu,
                  const std::vector<int> &load_type,
                  const std::vector<int> &load_node,
                  const std::vector<std::vector<std::complex<double>>> &u,
                  const std::vector<std::complex<double>> &i_ref) {
    for (auto &row : rhs) {
        std::fill(row.begin(), row.end(), std::complex<double>(0.0, 0.0));
    }

    for (size_t i = 0; i < load_type.size(); ++i) {
        int node_i = load_node[i];
        int type_i = load_type[i];
        if (type_i == CONST_POWER) {
            for (size_t j = 0; j < rhs.size(); ++j) {
                rhs[j][node_i] -= std::conj(load_pu[i] / u[j][node_i]);
            }
        } else if (type_i == CONST_CURRENT) {
            for (size_t j = 0; j < rhs.size(); ++j) {
                rhs[j][node_i] -= std::conj(
                    load_pu[i] * std::abs(u[j][node_i]) / u[j][node_i]);
            }
        } else if (type_i == CONST_IMPEDANCE) {
            for (size_t j = 0; j < rhs.size(); ++j) {
                rhs[j][node_i] -= std::conj(load_pu[i]) * u[j][node_i];
            }
        }
    }

    for (size_t j = 0; j < rhs.size(); ++j) {
        rhs[j].back() += i_ref[j];
    }
}

void solve_rhs_inplace_impl(
    const std::vector<int> &indptr_l, const std::vector<int> &indices_l,
    const std::vector<std::complex<double>> &data_l,
    const std::vector<int> &indptr_u, const std::vector<int> &indices_u,
    const std::vector<std::complex<double>> &data_u,
    std::vector<std::vector<std::complex<double>>> &rhs) {
    size_t size = rhs[0].size();
    // forward substitution
    for (size_t i = 0; i < size; ++i) {
        for (int index_j = indptr_l[i]; index_j < indptr_l[i + 1] - 1;
             ++index_j) {
            int j = indices_l[index_j];
            for (auto &row : rhs) {
                row[i] -= data_l[index_j] * row[j];
            }
        }
    }
    // backward substitution
    for (int i = size - 1; i >= 0; --i) {
        for (int index_j = indptr_u[i + 1] - 1; index_j > indptr_u[i];
             --index_j) {
            int j = indices_u[index_j];
            for (auto &row : rhs) {
                row[i] -= data_u[index_j] * row[j];
            }
        }
        int index_diag = indptr_u[i];
        for (auto &row : rhs) {
            row[i] /= data_u[index_diag];
        }
    }
}

double iterate_and_compare_impl(
    std::vector<std::vector<std::complex<double>>> &u,
    const std::vector<std::vector<std::complex<double>>> &rhs) {
    size_t size = u[0].size();
    double max_diff2 = 0.0;
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < u.size(); ++j) {
            std::complex<double> diff2 = rhs[j][i] - u[j][i];
            double diff2_val = std::norm(diff2);
            if (diff2_val > max_diff2) {
                max_diff2 = diff2_val;
            }
            u[j][i] = rhs[j][i];
        }
    }
    return max_diff2;
}

NAMESPACE_END
