#ifndef QUADRA_SPARSE_FACTORIZATION_HPP
#define QUADRA_SPARSE_FACTORIZATION_HPP
#pragma once

#include <iostream>
#include <stdexcept>
#include <string>

#include "../eigen/Eigen/Sparse"
#include "../eigen/Eigen/SparseCholesky"

namespace quadra {

/**
 * @brief Options controlling sparse matrix stabilization.
 *
 * Sparse mixed-effects models use the random-effect Hessian
 *
 * \f[
 *   H_{uu} = \frac{\partial^2 f(\theta, u)}{\partial u \partial u}.
 * \f]
 *
 * If factorization fails, Quadra can try
 *
 * \f[
 *   H_{uu}^{(\lambda)} = H_{uu} + \lambda I.
 * \f]
 */
struct SparseOptions {
    /// Initial diagonal jitter \f$\lambda\f$.
    double jitter_initial = 1e-12;

    /// Maximum number of adaptive jitter attempts.
    int jitter_max_attempts = 12;

    /// Print stabilization diagnostics.
    bool verbose = true;
};

/**
 * @brief Add diagonal jitter to a sparse matrix.
 *
 * This inserts diagonal elements when absent:
 *
 * \f[
 *   H_{ii} \leftarrow H_{ii} + \lambda.
 * \f]
 */
inline Eigen::SparseMatrix<double> add_diagonal_jitter(
    const Eigen::SparseMatrix<double>& H,
    double jitter)
{
    Eigen::SparseMatrix<double> H_reg = H;
    for (Eigen::Index i = 0; i < H_reg.rows(); ++i) {
        H_reg.coeffRef(i, i) += jitter;
    }
    H_reg.makeCompressed();
    return H_reg;
}

/**
 * @brief Factorize a sparse SPD matrix with adaptive jitter fallback.
 *
 * The unmodified matrix is tried first. Only if factorization fails does
 * Quadra try \f$H + \lambda I\f$.
 */
inline Eigen::SparseMatrix<double> factorize_with_adaptive_jitter(
    const Eigen::SparseMatrix<double>& H,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& solver,
    const char* context,
    const SparseOptions& options = SparseOptions())
{
    Eigen::SparseMatrix<double> H_factor = H;
    H_factor.makeCompressed();

    solver.compute(H_factor);
    if (solver.info() == Eigen::Success) {
        return H_factor;
    }

    double jitter = options.jitter_initial;
    for (int attempt = 0; attempt < options.jitter_max_attempts; ++attempt) {
        H_factor = add_diagonal_jitter(H, jitter);
        solver.compute(H_factor);

        if (solver.info() == Eigen::Success) {
            if (options.verbose) {
                std::cout << "Quadra: " << context
                          << " succeeded with diagonal jitter = "
                          << jitter << "\n";
            }
            return H_factor;
        }

        jitter *= 10.0;
    }

    throw std::runtime_error(std::string(context) + ": sparse factorization failed");
}

} // namespace quadra

#endif
