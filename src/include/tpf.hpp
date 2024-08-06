// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "./common.hpp"
#include "./utils.hpp"

NAMESPACE_BEGIN

// Minimal TPF impl

void set_load_pu(VectorComplex& load_pu, VectorScalar const& p_array, VectorScalar const& q_array) {
    std::transform(p_array.begin(), p_array.end(), q_array.begin(), load_pu.begin(),
                   [](Float p, Float q) { return Complex(p, q) / BASE_POWER; });
}

VectorComplex get_load_pu(VectorScalar const& p_specified, VectorScalar const& q_specified) {
    VectorComplex load_pu(p_specified.size());
    set_load_pu(load_pu, p_specified, q_specified);
    return load_pu;
}

void set_rhs(TensorComplex& rhs, VectorComplex const& load_pu, VectorInt const& load_type, VectorInt const& load_node,
             TensorComplex const& u, VectorComplex const& i_ref) {
    for (auto& row : rhs) {
        std::fill(row.begin(), row.end(), Complex(0.0, 0.0));
    }

    for (IDx i = 0; i < static_cast<IDx>(load_type.size()); ++i) {
        Int node_i = load_node[i];
        Int type_i = load_type[i];
        if (type_i == CONST_POWER) {
            for (IDx j = 0; j < static_cast<IDx>(rhs.size()); ++j) {
                rhs[j][node_i] -= std::conj(load_pu[i] / u[j][node_i]);
            }
        } else if (type_i == CONST_CURRENT) {
            for (IDx j = 0; j < static_cast<IDx>(rhs.size()); ++j) {
                rhs[j][node_i] -= std::conj(load_pu[i] * std::abs(u[j][node_i]) / u[j][node_i]);
            }
        } else if (type_i == CONST_IMPEDANCE) {
            for (IDx j = 0; j < static_cast<IDx>(rhs.size()); ++j) {
                rhs[j][node_i] -= std::conj(load_pu[i]) * u[j][node_i];
            }
        }
    }

    for (IDx j = 0; j < static_cast<IDx>(rhs.size()); ++j) {
        rhs[j].back() += i_ref[j];
    }
}

void solve_rhs_inplace(VectorInt const& indptr_l, VectorInt const& indices_l, VectorComplex const& data_l,
                       VectorInt const& indptr_u, VectorInt const& indices_u, VectorComplex const& data_u,
                       TensorComplex& rhs) {
    Int size = static_cast<IDx>(rhs[0].size());

    for (IDx i = 0; i < size; ++i) {
        for (IDx index_j = indptr_l[i]; index_j < indptr_l[i + 1] - 1; ++index_j) {
            Int j = indices_l[index_j];
            for (auto& row : rhs) {
                row[i] -= data_l[index_j] * row[j];
            }
        }
    }

    for (Int i = size - 1; i >= 0; --i) {
        for (Int index_j = indptr_u[i + 1] - 1; index_j > indptr_u[i]; --index_j) {
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

Float iterate_and_compare(TensorComplex& u, TensorComplex const& rhs) {
    Int size = static_cast<Int>(u[0].size());
    Float max_diff = 0.0;
    for (Int i = 0; i < size; ++i) {
        for (Int j = 0; j < static_cast<Int>(u.size()); ++j) {
            Complex diff = rhs[j][i] - u[j][i];
            Float diff_val = std::norm(diff);
            if (diff_val > max_diff) {
                max_diff = diff_val;
            }
            u[j][i] = rhs[j][i];
        }
    }
    return max_diff;
}

class TPF {
  public:
    TPF(PgmDataset const& input_data, Float const system_frequency, Int const n_node, Int const n_line)
        : _input_data(input_data), _n_node(n_node), _n_line(n_line), _system_frequency(system_frequency) {
        _node_org_to_reordered.resize(_n_node, -1);
    }

    ~TPF() = default;

    PgmResultType calculate_power_flow(PgmBatchDataset const& /*update_data*/, Int /*max_iteration*/ = 20,
                                       Float /*error_tolerance*/ = 1e-8) {
        return PgmResultType{};
    }

  private:
    void factorize_matrix() {
        SparsMatComplex y_matrix(_y_bus);
        y_matrix.makeCompressed();

        Eigen::SimplicialLLT<SparsMatComplex> solver;
        solver.analyzePattern(y_matrix);
        solver.factorize(y_matrix);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Matrix factorization failed!");
        }

        _l_matrix = solver.matrixL();
        _u_matrix = solver.matrixU();
    }

    void reorder_nodes(VectorInt const& reordered_node) {
        if (static_cast<Int>(reordered_node.size()) != _n_node) {
            throw std::invalid_argument("The graph is not connected!");
        }
        VectorInt reversed_node = reordered_node;
        std::reverse(reversed_node.begin(), reversed_node.end());
        _node_reordered_to_org = reversed_node;
        std::iota(_node_org_to_reordered.begin(), _node_org_to_reordered.end(), 0);
        for (size_t i = 0; i < reordered_node.size(); ++i) {
            _node_org_to_reordered[reordered_node[i]] = i;
        }
        _line_node_from = reorder_vector(_line_node_from);
        _line_node_to = reorder_vector(_line_node_to);
        _load_node = reorder_vector(_load_node);
        _source_node = _node_org_to_reordered[_source_node];
        assert(_source_node == _n_node - 1);
    }

    void graph_reorder() {} /* ToDo */

    void build_y_bus() {
        VectorComplex y_series(_n_line);
        for (Int i = 0; i < _n_line; ++i) {
            y_series[i] =
                1.0 / Complex(_input_data.at("line").at("r1").data(i), _input_data.at("line").at("x1").data(i));
        }
        VectorComplex y_shunt(_n_line);
        for (Int i = 0; i < _n_line; ++i) {
            y_shunt[i] = (0.5 * 2 * M_PI * _system_frequency) * _input_data.at("line").at("c1").data(i) *
                         Complex(0, 1 + _input_data.at("line").at("tan1").data(i));
        }
        VectorComplex all_y;
        all_y.insert(all_y.end(), y_series.begin(), y_series.end());
        all_y.insert(all_y.end(), y_shunt.begin(), y_shunt.end());
        all_y.insert(all_y.end(), y_shunt.begin(), y_shunt.end());

        SparsMatComplex y_branch(_n_line * 3, _n_line * 3);
        for (Int i = 0; i < _n_line * 3; ++i) {
            y_branch.insert(i, i) = all_y[i] / _y_base;
        }

        VectorInt ones(_n_line, 1);
        VectorInt incidence_entry;
        incidence_entry.insert(incidence_entry.end(), ones.begin(), ones.end());
        for (IDx i = 0; i < static_cast<IDx>(ones.size()); ++i) {
            incidence_entry.push_back(-1);
        }
        incidence_entry.insert(incidence_entry.end(), ones.begin(), ones.end());
        incidence_entry.insert(incidence_entry.end(), ones.begin(), ones.end());

        VectorInt incidence_i;
        incidence_i.insert(incidence_i.end(), _line_node_from.begin(), _line_node_from.end());
        incidence_i.insert(incidence_i.end(), _line_node_to.begin(), _line_node_to.end());
        incidence_i.insert(incidence_i.end(), _line_node_from.begin(), _line_node_from.end());
        incidence_i.insert(incidence_i.end(), _line_node_to.begin(), _line_node_to.end());

        VectorInt incidence_j;
        for (Int i = 0; i < _n_line; ++i) {
            incidence_j.push_back(i);
        }
        for (Int i = 0; i < _n_line; ++i) {
            incidence_j.push_back(i);
        }
        for (Int i = _n_line; i < _n_line * 2; ++i) {
            incidence_j.push_back(i);
        }
        for (Int i = _n_line * 2; i < _n_line * 3; ++i) {
            incidence_j.push_back(i);
        }

        SparsMatInt incidence_matrix_int(_n_node, _n_line * 3);
        for (size_t i = 0; i < incidence_entry.size(); ++i) {
            incidence_matrix_int.insert(incidence_i[i], incidence_j[i]) = incidence_entry[i];
        }

        SparsMatComplex incidence_matrix = incidence_matrix_int.cast<Complex>();
        _y_bus = (incidence_matrix * y_branch * incidence_matrix.transpose()).eval();
    }

    VectorInt reorder_vector(VectorInt const& vec) const {
        VectorInt reordered(vec.size());
        std::transform(vec.begin(), vec.end(), reordered.begin(),
                       [this](Int const node) { return _node_org_to_reordered[node]; });
        return reordered;
    }

    void pre_cache_calculation() {}

    PgmDataset _input_data;

    Int _n_node;
    Int _n_line;
    Int _source_node;

    Float _y_base;
    Float _system_frequency;

    VectorInt _node_reordered_to_org;
    VectorInt _node_org_to_reordered;
    VectorInt _line_node_from;
    VectorInt _line_node_to;
    VectorInt _load_node;
    VectorInt _load_type;

    SparsMatComplex _y_bus;
    SparsMatComplex _l_matrix;
    SparsMatComplex _u_matrix;
};

NAMESPACE_END
