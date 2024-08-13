// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "./common.hpp"
#include "./utils.hpp"

#include <omp.h>

NAMESPACE_BEGIN

// Minimal TPF impl

void set_load_pu(DenseMatComplex& load_pu, DenseMatReal const& p_matrix, DenseMatReal const& q_matrix) {
    assert(p_matrix.rows() == q_matrix.rows() && p_matrix.cols() == q_matrix.cols());
    load_pu.resize(p_matrix.rows(), p_matrix.cols());

    Float const inv_base_power = 1.0 / BASE_POWER;

#pragma omp parallel
    for (int j = 0; j < p_matrix.cols(); ++j) {
        for (int i = 0; i < p_matrix.rows(); ++i) {
            load_pu(i, j) = Complex(p_matrix(i, j), q_matrix(i, j)) * inv_base_power;
        }
    }
}

DenseMatComplex get_load_pu(PgmData const& update_data) {
    DenseMatReal const p_specified = update_data.at("p_specified").data;
    DenseMatReal const q_specified = update_data.at("q_specified").data;
    DenseMatComplex load_pu(p_specified.rows(), p_specified.cols());

    set_load_pu(load_pu, p_specified, q_specified);
    return load_pu;
}

void set_rhs(SparseMatComplex& rhs, DenseMatComplex const& load_pu, VectorInt const& load_type,
             VectorInt const& load_node, SparseMatComplex const& u, Complex const& i_ref) {
    DenseMatComplex u_dense = u;
    DenseMatComplex rhs_dense = rhs;
    rhs_dense.setZero();

    for (IDx i = 0; i < static_cast<IDx>(load_type.size()); ++i) {
        Int node_i = load_node[i];
        Int type_i = load_type[i];
        switch (type_i) {
        case CONST_POWER:
            for (IDx j = 0; j < static_cast<IDx>(rhs.innerSize()); ++j) {
                rhs_dense(j, node_i) -= std::conj(load_pu(j, i) / u_dense(j, node_i));
            }
            break;

        case CONST_CURRENT:
            for (IDx j = 0; j < static_cast<IDx>(rhs.innerSize()); ++j) {
                rhs_dense(j, node_i) -= std::conj(load_pu(j, i) * std::abs(u_dense(j, node_i)) / u_dense(j, node_i));
            }
            break;

        case CONST_IMPEDANCE:
            for (IDx j = 0; j < static_cast<IDx>(rhs.innerSize()); ++j) {
                rhs_dense(j, node_i) -= std::conj(load_pu(j, i)) * u_dense(j, node_i);
            }
            break;

        default:
            throw std::runtime_error("Unexpected type_i value: " + std::to_string(type_i));
            break;
        }
    }

    for (IDx j = 0; j < rhs_dense.rows(); ++j) {
        rhs_dense(j, rhs_dense.cols() - 1) += i_ref;
    }
    rhs = rhs_dense.sparseView();
}

class TPF {
  public:
    TPF(PgmData const& input_data, Float const system_frequency)
        : _input_data(input_data), _system_frequency(system_frequency), _solver() {
        _n_node = static_cast<Int>(input_data.at("node").data.rows());
        _n_line = static_cast<Int>(input_data.at("line").data.rows());
        _node_org_to_reordered.resize(_n_node, -1);

        _u_rated = input_data.at("node").data.col(1)[0];
        assert(_u_rated != 0.0);
        _y_base = BASE_POWER / (_u_rated * _u_rated);

        auto retrieve_from_data = [](PgmData const& in_data, char const* component, Int col) {
            auto const& data_col = in_data.at(component).data.col(col);
            return VectorInt(data_col.data(), data_col.data() + data_col.size());
        };

        _line_node_from = retrieve_from_data(input_data, "line", 1);
        _line_node_to = retrieve_from_data(input_data, "line", 2);
        _load_node = retrieve_from_data(input_data, "sym_load", 1);
        _load_type = retrieve_from_data(input_data, "sym_load", 3);

        _source_node = static_cast<Int>(input_data.at("source").data(0, 1));

        pre_cache_calculation();
    }

    ~TPF() = default;

    void printDenseMat(DenseMatReal const& mat, std::string const& name) {
        std::cout << name << ": \n[ " << std::endl;
        for (int i = 0; i < mat.rows(); ++i) {
            std::cout << "[ ";
            for (int j = 0; j < mat.cols(); ++j) {
                std::cout << mat(i, j) << " ";
            }
            std::cout << " ]" << std::endl;
        }
        std::cout << " ]" << std::endl;
    }

    PgmResultType calculate_power_flow(PgmData const& update_data, Int max_iteration = 20,
                                       Float error_tolerance = 1e-8) {
        pre_cache_calculation();

        Int const n_steps = static_cast<Int>(update_data.at("p_specified").data.rows());

        // Initialize load_pu, u, and rhs
        DenseMatComplex const load_pu = get_load_pu(update_data);
        DenseMatComplex u_dense(n_steps, _n_node);
        u_dense.setConstant(_u_ref);
        SparseMatComplex u = u_dense.sparseView();
        SparseMatComplex rhs(n_steps, _n_node);

        // Iterate
        for (Int iter = 0; iter < max_iteration; ++iter) {
            set_rhs(rhs, load_pu, _load_type, _load_node, u, _i_ref);
            solve_rhs_inplace(rhs);
            Float max_diff = iterate_and_compare(u, rhs);

            // Early out, 'converged'
            if (max_diff < error_tolerance * error_tolerance) {
                break;
            }

            if (iter == max_iteration - 1) {
                throw std::runtime_error("The power flow calculation does not converge! Max diff: " +
                                         std::to_string(max_diff));
            }
        }

        // Reorder back to original
        DenseMatReal u_pu(n_steps, _n_node);
        DenseMatReal u_angle(n_steps, _n_node);
#pragma omp parallel
        for (Int i = 0; i < n_steps; ++i) {
            for (Int j = 0; j < _n_node; ++j) {
                Complex u_value = u.coeff(i, _node_org_to_reordered[j]);
                u_pu.coeffRef(i, j) = std::abs(u_value);
                u_angle.coeffRef(i, j) = std::arg(u_pu.coeff(i, j));
            }
        }

        return produce_result(u_pu, u_angle);
    }

  private:
    // Matrix operations inside iterations ========================================
    void solve_rhs_inplace(SparseMatComplex& rhs) const {
        DenseMatComplex rhs_dense = rhs;
        DenseMatComplex rhs_dense_t = rhs_dense.transpose();

        _solver.matrixL().solveInPlace(rhs_dense_t);
        _solver.matrixU().solveInPlace(rhs_dense_t);

        rhs_dense = rhs_dense_t.transpose();
        rhs = rhs_dense.sparseView();
    }

    Float iterate_and_compare(SparseMatComplex& u, SparseMatComplex const& rhs) {
        DenseMatComplex u_dense = u;
        DenseMatComplex rhs_dense = rhs;

        Int size = static_cast<Int>(u_dense.cols());
        Float max_diff = 0.0;

#pragma omp parallel
        for (Int i = 0; i < size; ++i) {
            for (Int j = 0; j < static_cast<Int>(u_dense.rows()); ++j) {
                std::complex<double> diff = rhs_dense(j, i) - u_dense(j, i);
                Float diff_val = std::norm(diff); // should have been squared norm
                if (diff_val > max_diff) {
                    max_diff = diff_val;
                }
                u_dense(j, i) = rhs_dense(j, i);
            }
        }

        u = u_dense.sparseView();

        return max_diff;
    }

    // graph ordering ==============================================================
    void reorder_nodes(VectorInt const& reordered_node) {
        auto const& reorder_nodes_size = static_cast<Int>(reordered_node.size());
        if (reorder_nodes_size != _n_node) {
            throw std::invalid_argument("The graph is not connected!");
        }
        VectorInt reversed_node = reordered_node;
        std::reverse(reversed_node.begin(), reversed_node.end());
        _node_reordered_to_org = reversed_node;
        std::iota(_node_org_to_reordered.begin(), _node_org_to_reordered.end(), 0);
        for (Int i = 0; i < reorder_nodes_size; ++i) {
            _node_org_to_reordered[reordered_node[i]] = reorder_nodes_size - 1 - i;
        }
        _line_node_from = reorder_vector(_line_node_from);
        _line_node_to = reorder_vector(_line_node_to);
        _load_node = reorder_vector(_load_node);
        _source_node = _node_org_to_reordered[_source_node];
        assert(_source_node == _n_node - 1);
    }

    VectorInt reorder_vector(VectorInt const& vec) {
        VectorInt reordered_vec(vec.size());
        for (Int i = 0; i < static_cast<Int>(vec.size()); ++i) {
            reordered_vec[i] = _node_org_to_reordered[vec[i]];
        }
        return reordered_vec;
    }

    struct DfsOrderRecorder : public boost::default_dfs_visitor {
        DfsOrderRecorder(VectorInt& order) : order(order) {}

        template <typename Vertex, typename Graph> void discover_vertex(Vertex u, const Graph& /*g*/) const {
            order.push_back(u);
        }

        VectorInt& order;
    };

    VectorInt depth_first_order(SparseMatInt const& connection_array, Int start_node) {
        using namespace boost;

        typedef adjacency_list<vecS, vecS, undirectedS> Graph;

        Graph g;
        for (Int k = 0; k < static_cast<Int>(connection_array.outerSize()); ++k) {
            for (SparseMatInt::InnerIterator it(connection_array, k); it; ++it) {
                add_edge(it.row(), it.col(), g);
            }
        }

        VectorInt dfs_order;
        DfsOrderRecorder vis(dfs_order);
        std::vector<default_color_type> color_map(num_vertices(g));
        depth_first_search(g, visitor(vis).root_vertex(start_node).color_map(&color_map[0]));

        return dfs_order;
    }

    void graph_reorder() {
        VectorInt edge_i(_line_node_from);
        edge_i.insert(edge_i.end(), _line_node_to.begin(), _line_node_to.end());

        VectorInt edge_j(_line_node_to);
        edge_j.insert(edge_j.end(), _line_node_from.begin(), _line_node_from.end());

        SparseMatInt connection_array(_n_node, _n_node);
        std::vector<Eigen::Triplet<Int>> tripletList;
        tripletList.reserve(edge_i.size());
        for (Int k = 0; k < static_cast<Int>(edge_i.size()); ++k) {
            tripletList.emplace_back(edge_i[k], edge_j[k], Int(1));
        }
        connection_array.setFromTriplets(tripletList.begin(), tripletList.end());

        VectorInt reordered_node = depth_first_order(connection_array, _source_node);

        reorder_nodes(reordered_node);
    }

    // matrix factorizatio and construction =======================================
    void factorize_matrix() {
        SparseMatComplex y_matrix(_y_bus);
        y_matrix.makeCompressed();

        _solver.analyzePattern(y_matrix);
        _solver.factorize(y_matrix);

        if (_solver.info() != Eigen::Success) {
            throw std::runtime_error("Matrix factorization failed!");
        }
        _is_initialized = true;
    }

    void build_y_bus() {
        VectorComplex y_series(_n_line);
        for (Int i = 0; i < _n_line; ++i) {
            Float r1 = _input_data.at("line").data(i, COL_R1);
            Float x1 = _input_data.at("line").data(i, COL_X1);
            y_series[i] = 1.0 / Complex(r1, x1);
        }
        VectorComplex y_shunt(_n_line);
        for (Int i = 0; i < _n_line; ++i) {
            Float c1 = _input_data.at("line").data(i, COL_C1);
            Float tan1 = _input_data.at("line").data(i, COL_TAN1);
            y_shunt[i] = (0.5 * 2 * pi * _system_frequency) * c1 * Complex(0, 1 + tan1);
        }
        VectorComplex all_y;
        all_y.insert(all_y.end(), y_series.begin(), y_series.end());
        all_y.insert(all_y.end(), y_shunt.begin(), y_shunt.end());
        all_y.insert(all_y.end(), y_shunt.begin(), y_shunt.end());

        SparseMatComplex y_branch(_n_line * 3, _n_line * 3);
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

        SparseMatInt incidence_matrix_int(_n_node, _n_line * 3);
        for (size_t i = 0; i < incidence_entry.size(); ++i) {
            incidence_matrix_int.insert(incidence_i[i], incidence_j[i]) = incidence_entry[i];
        }

        SparseMatComplex incidence_matrix = incidence_matrix_int.cast<Complex>();
        _y_bus = (incidence_matrix * y_branch * incidence_matrix.transpose()).eval();
    }

    // initialization =============================================================
    void pre_cache_calculation() {
        if (_node_reordered_to_org.empty()) {
            graph_reorder();
        }
        if (_y_bus.nonZeros() == 0) {
            build_y_bus();
        }
        if (!_is_initialized) {
            factorize_matrix();
        }
        // u variable, flat start as u_ref
        auto const& source_node_data = _input_data.at("source").data;
        _u_ref = Complex(source_node_data(0, 0), 0.0);
        _i_ref = _y_base * _u_ref;
    }

    PgmResultType produce_result(DenseMatReal const& u_pu, DenseMatReal const& u_angle) {
        PgmData u_pu_data;
        u_pu_data["value"].data = u_pu;

        PgmData u_angle_data;
        u_angle_data["value"].data = u_angle;

        PgmDataset node_result;
        node_result["u_pu"] = u_pu_data;
        node_result["u_angle"] = u_angle_data;

        PgmResultType result;
        result["node"] = node_result;

        return result;
    }

    PgmData _input_data;

    Int _n_node;
    Int _n_line;
    Int _source_node;

    bool _is_initialized = false;

    Float _u_rated;
    Float _y_base;
    Float _system_frequency;

    Complex _u_ref;
    Complex _i_ref;

    VectorInt _node_reordered_to_org;
    VectorInt _node_org_to_reordered;
    VectorInt _line_node_from;
    VectorInt _line_node_to;
    VectorInt _load_node;
    VectorInt _load_type;

    SparseMatComplex _y_bus;
    Eigen::SparseLU<SparseMatComplex> _solver;
};

NAMESPACE_END
