#ifndef QUADRA_LAPLACE_HPP
#define QUADRA_LAPLACE_HPP
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include <unordered_map>
#include <set>
#include <iomanip>
#include "../model/parameter.hpp"
#include "autodiff.hpp"
#include "evaluation.hpp"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/SparseCholesky"

namespace quadra
{

    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    //==============================
    // Build fixed index map
    //==============================
    inline std::vector<int>
    build_fixed_index(const ParameterVector &params)
    {
        std::vector<int> idx;
        for (size_t i = 0; i < params.params.size(); ++i)
        {
            if (!params.params[i].is_random)
                idx.push_back(i);
        }
        return idx;
    }

    //==============================
    // Build random index map
    //==============================
    inline std::vector<int>
    build_random_index(const ParameterVector &params)
    {
        std::vector<int> idx;
        for (size_t i = 0; i < params.params.size(); ++i)
        {
            if (params.params[i].is_random)
                idx.push_back(i);
        }
        return idx;
    }

    std::vector<double> build_u_init_from_cache(
        const std::vector<int> &random_idx)
    {
        return std::vector<double>(random_idx.size(), 0.0);
    }

    inline void inject_fixed_params(
        const Eigen::VectorXd &x,
        ParameterVector &params,
        const std::vector<int> &fixed_idx)
    {
        assert(x.size() == fixed_idx.size());

        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            const int idx = fixed_idx[k];
            params.params[idx].value = x[k];
        }
    }

    template <typename Scalar>
    inline void inject_fixed_params(
        const std::vector<Scalar> &x_ad,
        std::vector<Scalar> &p,
        const std::vector<int> &fixed_idx)
    {
        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            p[fixed_idx[k]] = x_ad[k];
        }
    }

    template <typename Scalar>
    inline void inject_fixed_params(
        const Eigen::VectorXd &x,
        std::vector<Scalar> &p,
        const std::vector<int> &fixed_idx)
    {
        for (size_t k = 0; k < fixed_idx.size(); ++k)
            p[fixed_idx[k]] = Scalar(x[k]);
    }

    inline void inject_random_params(
        const std::vector<double> &u,
        ParameterVector &params,
        const std::vector<int> &random_idx)
    {
        assert(u.size() == random_idx.size());

        for (size_t k = 0; k < random_idx.size(); ++k)
        {
            const int idx = random_idx[k];
            params.params[idx].value = u[k];
        }
    }

    template <typename Scalar>
    inline void inject_random_params(
        const std::vector<Scalar> &u_ad,
        std::vector<Scalar> &p,
        const std::vector<int> &random_idx)
    {
        for (size_t k = 0; k < random_idx.size(); ++k)
        {
            p[random_idx[k]] = u_ad[k];
        }
    }

    template <typename Scalar>
    inline void inject_random_params(
        const std::vector<double> &u,
        std::vector<Scalar> &p,
        const std::vector<int> &random_idx)
    {
        for (size_t k = 0; k < random_idx.size(); ++k)
        {
            p[random_idx[k]] = Scalar(u[k]);
        }
    }

    template <typename T>
    std::vector<T> pack_params(
        const std::vector<T> &u,
        const std::vector<T> &x,
        const ParameterVector &params,
        const std::vector<int> &random_idx,
        const std::vector<int> &fixed_idx)
    {
        std::vector<T> p(params.params.size());

        int u_k = 0;
        int x_k = 0;

        for (size_t i = 0; i < p.size(); i++)
        {
            if (params.params[i].is_random)
            {
                p[i] = u[u_k++];
            }
            else
            {
                p[i] = x[x_k++];
            }
        }

        return p;
    }

    //==================================================
    // Laplace-local Hessian pattern representation
    //==================================================
    // Do not name this HessianPattern. autodiff.hpp may define a
    // graph-level HessianPattern helper for ADGraph sparsity discovery.
    // Keeping the Laplace cache as SparseHessianPattern avoids redefinition
    // errors and keeps this file independent of the exact autodiff helper API.
    using SparseHessianPattern = std::vector<std::pair<int, int>>;

    inline std::unordered_map<int, SparseHessianPattern> &laplace_pattern_cache()
    {
        static std::unordered_map<int, SparseHessianPattern> cache;
        return cache;
    }

    //==================================================
    // Discover Hessian sparsity from had::ADGraph
    //==================================================
    // This replaces the older dense pattern probe. It reads the sparse
    // edge-pushed Hessian storage that had::PropagateAdjoint() has already
    // populated inside scope.backward(nll).
    //
    // NOTE: this is still a numeric sparsity pattern. If a structurally
    // nonzero Hessian entry evaluates to exactly zero at the discovery point,
    // it can be missed. Diagonals are included by default for Newton stability.
    inline SparseHessianPattern discover_pattern_from_graph(
        const std::vector<AD> &p_full,
        const std::vector<int> &random_idx,
        bool symmetric = true,
        bool include_diagonal = true,
        double tol = 1e-12)
    {
        std::cout << "Quadra: Discovering Hessian pattern from AD graph for "
                  << random_idx.size() << " random variables ...\n";

        const int n = static_cast<int>(random_idx.size());
        SparseHessianPattern pattern;

        if (n == 0 || had::g_ADGraph == nullptr)
            return pattern;

        std::unordered_map<had::VertexId, int> random_var_to_local;
        random_var_to_local.reserve(static_cast<size_t>(n));

        for (int local = 0; local < n; ++local)
        {
            const int full_index = random_idx[static_cast<size_t>(local)];
            random_var_to_local.emplace(p_full[full_index].varId, local);
        }

        std::set<std::pair<int, int>> unique_pairs;

        if (include_diagonal)
        {
            for (int i = 0; i < n; ++i)
                unique_pairs.emplace(i, i);
        }
        else
        {
            for (int i = 0; i < n; ++i)
            {
                const int full_index = random_idx[static_cast<size_t>(i)];
                const had::VertexId vi = p_full[full_index].varId;

                if (vi < had::g_ADGraph->selfSoEdges.size() &&
                    std::abs(had::g_ADGraph->selfSoEdges[vi]) > tol)
                {
                    unique_pairs.emplace(i, i);
                }
            }
        }

        // had stores an off-diagonal Hessian entry in soEdges[max_id]
        // under key min_id. Walk the graph-level sparse storage and retain
        // only entries where both endpoints are random-effect variables.
        for (had::VertexId hi = 0;
             hi < static_cast<had::VertexId>(had::g_ADGraph->soEdges.size());
             ++hi)
        {
            auto hi_it = random_var_to_local.find(hi);
            if (hi_it == random_var_to_local.end())
                continue;

            const int i = hi_it->second;
            const auto &tree = had::g_ADGraph->soEdges[hi];

            for (const auto &node : tree.nodes)
            {
                if (std::abs(node.val) <= tol)
                    continue;

                auto lo_it = random_var_to_local.find(node.key);
                if (lo_it == random_var_to_local.end())
                    continue;

                const int j = lo_it->second;
                unique_pairs.emplace(i, j);
                if (symmetric)
                    unique_pairs.emplace(j, i);
            }
        }

        pattern.reserve(unique_pairs.size());
        for (const auto &ij : unique_pairs)
            pattern.emplace_back(ij.first, ij.second);

        std::cout << "Quadra: Model structure aware now => Hessian pattern has " << pattern.size() << " entries.\n";
        return pattern;
    }

    inline const SparseHessianPattern &get_pattern(
        const ADScope &,
        const std::vector<AD> &p_full,
        const std::vector<int> &random_idx)
    {
        const int n = static_cast<int>(random_idx.size());
        auto &cache = laplace_pattern_cache();

        auto it = cache.find(n);
        if (it != cache.end())
            return it->second;

        auto pattern = discover_pattern_from_graph(p_full, random_idx);
        auto res = cache.emplace(n, std::move(pattern));
        return res.first->second;
    }

    inline SparseHessianPattern banded_hessian_pattern(int n, int bandwidth)
    {
        SparseHessianPattern pattern;

        for (int i = 0; i < n; ++i)
        {
            int j0 = std::max(0, i - bandwidth);
            int j1 = std::min(n - 1, i + bandwidth);

            for (int j = j0; j <= j1; ++j)
                pattern.emplace_back(i, j);
        }

        return pattern;
    }

    inline SparseHessianPattern dense_hessian_pattern(int n)
    {
        SparseHessianPattern pattern;
        pattern.reserve(n * n);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                pattern.emplace_back(i, j);

        return pattern;
    }

    inline Eigen::SparseMatrix<double> extract_sparse_hessian(
        const ADScope &scope,
        const std::vector<AD> &p_full,
        const std::vector<int> &random_idx,
        const SparseHessianPattern &pattern,
        double drop_tol = 1e-12)
    {
        const int n = (int)random_idx.size();

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(pattern.size());

        for (const auto &[i, j] : pattern)
        {
            double hij = scope.hess(
                p_full[random_idx[i]],
                p_full[random_idx[j]]);
            if (std::abs(hij) > drop_tol)
                triplets.emplace_back(i, j, hij);
        }

        Eigen::SparseMatrix<double> H(n, n);
        H.setFromTriplets(triplets.begin(), triplets.end());
        return H;
    }

    inline double sparse_logdet_ldlt(
        const Eigen::SparseMatrix<double> &H)
    {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(H);

        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Sparse LDLT factorization failed");
        }

        const auto &D = solver.vectorD();

        double logdet = 0.0;

        for (int i = 0; i < D.size(); ++i)
        {
            if (D[i] <= 0.0)
            {
                throw std::runtime_error("Sparse Hessian is not positive definite");
            }

            logdet += std::log(D[i]);
        }

        return logdet;
    }

    inline double sparse_logdet_llt(
        const Eigen::SparseMatrix<double> &H)
    {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(H);

        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Sparse LLT factorization failed");
        }

        Eigen::SparseMatrix<double> L = solver.matrixL();

        double logdet = 0.0;

        for (int k = 0; k < L.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
            {
                if (it.row() == it.col())
                {
                    if (it.value() <= 0.0)
                    {
                        throw std::runtime_error("Non-positive Cholesky diagonal");
                    }

                    logdet += 2.0 * std::log(it.value());
                }
            }
        }

        return logdet;
    }
    //==================================================
    // Solve for random effects u* via Newton
    //==================================================
    template <typename Model>
    std::vector<double> solve_random_effects_laplace(
        Model &model,
        ParameterVector &params,
        const Eigen::VectorXd &x,
        const std::vector<int> &fixed_idx,
        const std::vector<int> &random_idx,
        had::ADGraph &graph)
    {
        const int max_iter = 20;
        const double tol = 1e-8;

        std::vector<double> u(random_idx.size(), 0.0);

        for (int iter = 0; iter < max_iter; ++iter)
        {

            // --------------------------------------------------
            // AD scope: binds graph, clears graph
            // --------------------------------------------------
            ADScope scope(graph);

            // --------------------------------------------------
            // Build full AD parameter vector
            // --------------------------------------------------
            std::vector<AD> p_full;
            p_full.reserve(params.size());

            for (int i = 0; i < params.size(); ++i)
            {
                p_full.emplace_back(AD(0.0));
            }

            inject_fixed_params(x, p_full, fixed_idx);
            inject_random_params(u, p_full, random_idx);

            // --------------------------------------------------
            // Forward pass
            // --------------------------------------------------
            AD nll = model(p_full);

            // --------------------------------------------------
            // Reverse pass
            // --------------------------------------------------
            scope.backward(nll);

            // --------------------------------------------------
            // Gradient wrt random effects
            // --------------------------------------------------
            Eigen::VectorXd g(random_idx.size());

            for (size_t i = 0; i < random_idx.size(); ++i)
            {
                g[i] = scope.grad(p_full[random_idx[i]]);
            }

                        if (g.norm() < tol)
            {
                return u;
            }

            // --------------------------------------------------
            // Sparse Hessian wrt random effects
            // --------------------------------------------------
            const auto &pattern = get_pattern(scope, p_full, random_idx);

            Eigen::SparseMatrix<double> H =
                extract_sparse_hessian(scope, p_full, random_idx, pattern);

            // Optional diagnostics
            // std::cout << "pattern nnz = " << H.nonZeros() << "\n";

            // --------------------------------------------------
            // Sparse Newton solve: H step = g
            // --------------------------------------------------
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(H);

            if (solver.info() != Eigen::Success)
            {
                throw std::runtime_error(
                    "Sparse Hessian factorization failed in solve_random_effects_laplace");
            }

            Eigen::VectorXd step = solver.solve(g);

            if (solver.info() != Eigen::Success)
            {
                throw std::runtime_error(
                    "Sparse Hessian solve failed in solve_random_effects_laplace");
            }
            std::cout << "Newton: "
                      << "inner iter = " << std::setw(3) << iter
                      << ", fx = " << std::setw(14) << std::fixed << std::setprecision(6) << nll.val
                      << ", |grad| = " << std::setw(12) << std::fixed << std::setprecision(6) << g.norm()
                      << "\n";
            // --------------------------------------------------
            // Newton update
            // --------------------------------------------------
            for (size_t i = 0; i < u.size(); ++i)
            {
                u[i] -= step[i];
            }
        }

        return u;
    }

    template <typename Model>
    struct LaplaceResult
    {
        double value;
        std::vector<double> grad_x;
        std::vector<double> grad_u;
    };

    template <typename Model>
    LaplaceResult<Model> laplace_eval_at_u_star(
        Model &model,
        ParameterVector &params,
        const std::vector<int> &fixed_idx,
        const std::vector<int> &random_idx,
        const Eigen::VectorXd &x,
        const std::vector<double> &u_star,
        had::ADGraph &graph)
    {
        ADScope scope(graph);

        using Result = LaplaceResult<Model>;
        Result res;

        std::vector<AD> p_full;
        p_full.reserve(params.size());

        for (int i = 0; i < params.size(); ++i)
        {
            p_full.emplace_back(AD(0.0));
        }

        inject_fixed_params(x, p_full, fixed_idx);
        inject_random_params(u_star, p_full, random_idx);

        AD nll = model(p_full);

        scope.backward(nll);

        res.grad_x.resize(fixed_idx.size());
        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            res.grad_x[k] = scope.grad(p_full[fixed_idx[k]]);
        }

        res.grad_u.resize(random_idx.size());
        for (size_t k = 0; k < random_idx.size(); ++k)
        {
            res.grad_u[k] = scope.grad(p_full[random_idx[k]]);
        }

        const auto &pattern = get_pattern(scope, p_full, random_idx);

        Eigen::SparseMatrix<double> H =
            extract_sparse_hessian(scope, p_full, random_idx, pattern);

        double logdet = sparse_logdet_ldlt(H);
        // Or, if vectorD() is unavailable:
        // double logdet = sparse_logdet_llt(H);

        res.value = value_of(nll) + 0.5 * logdet;

        return res;
    }

#ifndef QUADRA_USE_ORIGINAL_HAD
    //==================================================
    // Optional third-order directional diagnostic.
    // This evaluates D^k f(x)[direction,...] for k = 0,1,2,3
    // using the scalar-templated model path. It is intentionally
    // separate from LBFGS/Laplace so it can be enabled only when needed.
    //==================================================
    template <typename Model>
    ThirdDirectionalResult third_directional_fixed_effects(
        Model &model,
        const Eigen::VectorXd &x,
        const Eigen::VectorXd &direction)
    {
        if (x.size() != direction.size())
            throw std::invalid_argument("third_directional_fixed_effects: x and direction must have same size");

        std::vector<double> xv(static_cast<size_t>(x.size()));
        std::vector<double> dv(static_cast<size_t>(direction.size()));
        for (int i = 0; i < x.size(); ++i)
        {
            xv[static_cast<size_t>(i)] = x[i];
            dv[static_cast<size_t>(i)] = direction[i];
        }

        return evaluate_third_directional(
            [&](const std::vector<AD3> &x_ad3) -> AD3
            {
                return model(x_ad3);
            },
            xv,
            dv);
    }
#endif

} // namespace quadra

#endif // QUADRA_LAPLACE_HPP
