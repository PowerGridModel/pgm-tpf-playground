#pragma once

#include "./common.hpp"
#include "./utils.hpp"

NAMESPACE_BEGIN

// Minimal TPF impl ==============================

void set_load_pu(VectorComplex& load_pu, VectorScalar const& p_array,
                 VectorScalar const& q_array) {
    std::transform(p_array.begin(), p_array.end(), q_array.begin(),
                   load_pu.begin(), [](Float p, Float q) {
                       return Complex(p, q) / BASE_POWER;
                   });
}

VectorComplex get_load_pu(VectorScalar const& p_specified,
                          VectorScalar const& q_specified) {
    VectorComplex load_pu(p_specified.size());
    set_load_pu(load_pu, p_specified, q_specified);
    return load_pu;
}

void set_rhs_impl(TensorComplex& rhs, VectorComplex const& load_pu,
                  VectorInt const& load_type, VectorInt const& load_node,
                  TensorComplex const& u, VectorComplex const& i_ref) {
    for (auto& row : rhs) {
        std::fill(row.begin(), row.end(), Complex(0.0, 0.0));
    }

    for (Int i = 0; i < static_cast<Int>(load_type.size()); ++i) {
        Int node_i = load_node[i];
        Int type_i = load_type[i];
        if (type_i == CONST_POWER) {
            for (Int j = 0; j < static_cast<Int>(rhs.size()); ++j) {
                rhs[j][node_i] -= std::conj(load_pu[i] / u[j][node_i]);
            }
        } else if (type_i == CONST_CURRENT) {
            for (Int j = 0; j < static_cast<Int>(rhs.size()); ++j) {
                rhs[j][node_i] -= std::conj(
                    load_pu[i] * std::abs(u[j][node_i]) / u[j][node_i]);
            }
        } else if (type_i == CONST_IMPEDANCE) {
            for (Int j = 0; j < static_cast<Int>(rhs.size()); ++j) {
                rhs[j][node_i] -= std::conj(load_pu[i]) * u[j][node_i];
            }
        }
    }

    for (Int j = 0; j < static_cast<Int>(rhs.size(); ++j) {
        rhs[j].back() += i_ref[j];
    }
}

void solve_rhs_inplace_impl(VectorInt const& indptr_l,
                            VectorInt const& indices_l,
                            VectorComplex const& data_l,
                            VectorInt const& indptr_u,
                            VectorInt const& indices_u, VectorComplex& data_u,
                            TensorComplex& rhs) {
    Int size = static_cast<Int>(rhs[0].size());
    // forward substitution
    for (Int i = 0; i < size; ++i) {
        for (Int index_j = indptr_l[i]; index_j < indptr_l[i + 1] - 1;
             ++index_j) {
            Int j = indices_l[index_j];
            for (auto& row : rhs) {
                row[i] -= data_l[index_j] * row[j];
            }
        }
    }
    // backward substitution
    for (Int i = size - 1; i >= 0; --i) {
        for (Int index_j = indptr_u[i + 1] - 1; index_j > indptr_u[i];
             --index_j) {
            Int j = indices_u[index_j];
            for (auto& row : rhs) {
                row[i] -= data_u[index_j] * row[j];
            }
        }
        Int index_diag = indptr_u[i];
        for (auto& row : rhs) {
            row[i] /= data_u[index_diag];
        }
    }
}

Float iterate_and_compare_impl(TensorComplex& u, TensorComplex const& rhs) {
    Int size = static_cast<Int>(u[0].size());
    Float max_diff2 = 0.0;
    for (Int i = 0; i < size; ++i) {
        for (Int j = 0; j < static_cast<Int>(u.size()); ++j) {
            Complex diff2 = rhs[j][i] - u[j][i];
            Float diff2_val = std::norm(diff2);
            if (diff2_val > max_diff2) {
                max_diff2 = diff2_val;
            }
            u[j][i] = rhs[j][i];
        }
    }
    return max_diff2;
}

NAMESPACE_END
