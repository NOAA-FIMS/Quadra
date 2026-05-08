#ifndef QUADRA_SPARSE_TRACE_HPP
#define QUADRA_SPARSE_TRACE_HPP
#pragma once

#include <random>
#include <stdexcept>

#include "../eigen/Eigen/Dense"
#include "../eigen/Eigen/Sparse"

namespace quadra {

/**
 * @brief Options for evaluating trace terms.
 *
 * Laplace gradients require
 *
 * \f[
 *   \operatorname{tr}\left(H^{-1}\dot H\right),
 * \f]
 *
 * where \f$H = H_{uu}\f$ and \f$\dot H = D H_{uu}[v]\f$.
 */
struct TraceOptions {
    /// Use Hutchinson stochastic trace estimation.
    bool use_hutchinson_trace = true;

    /// Number of Rademacher probe vectors.
    int hutchinson_probes = 8;

    /// RNG seed for reproducible stochastic traces.
    unsigned int hutchinson_seed = 12345;
};

/**
 * @brief Compute or estimate \f$\operatorname{tr}(H^{-1}\dot H)\f$.
 *
 * Dense exact trace:
 *
 * \f[
 *   \operatorname{tr}(H^{-1}\dot H).
 * \f]
 *
 * Hutchinson trace estimator:
 *
 * \f[
 *   \operatorname{tr}(A)
 *   =
 *   \mathbb{E}_z[z^\top A z],
 *   \quad z_i \in \{-1, 1\}.
 * \f]
 */
template <typename SolverType>
double trace_hinv_hdot(
    SolverType& solver,
    const Eigen::SparseMatrix<double>& Hdot,
    const TraceOptions& options = TraceOptions())
{
    if (Hdot.rows() != Hdot.cols()) {
        throw std::invalid_argument("trace_hinv_hdot: Hdot must be square");
    }

    const Eigen::Index n = Hdot.rows();

    if (!options.use_hutchinson_trace) {
        Eigen::MatrixXd rhs = Eigen::MatrixXd(Hdot);
        Eigen::MatrixXd X = solver.solve(rhs);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("trace_hinv_hdot: dense trace solve failed");
        }

        return X.diagonal().sum();
    }

    std::mt19937 rng(options.hutchinson_seed);
    std::uniform_int_distribution<int> rademacher(0, 1);

    double trace_est = 0.0;

    for (int sample = 0; sample < options.hutchinson_probes; ++sample) {
        Eigen::VectorXd z(n);
        for (Eigen::Index i = 0; i < n; ++i) {
            z[i] = (rademacher(rng) == 0) ? -1.0 : 1.0;
        }

        Eigen::VectorXd y = Hdot * z;
        Eigen::VectorXd x = solver.solve(y);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("trace_hinv_hdot: Hutchinson sparse solve failed");
        }

        trace_est += z.dot(x);
    }

    return trace_est / static_cast<double>(options.hutchinson_probes);
}

} // namespace quadra

#endif
