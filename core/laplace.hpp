#ifndef QUADRA_LAPLACE_HPP
#define QUADRA_LAPLACE_HPP
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "../model/parameter.hpp"
#include "autodiff.hpp"
#include "evaluation.hpp"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/SparseCholesky"

DECLARE_ADGRAPH();

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

    inline void inject_fixed_params(
        const std::vector<had::AReal> &x_ad,
        std::vector<had::AReal> &p,
        const std::vector<int> &fixed_idx)
    {
        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            p[fixed_idx[k]] = x_ad[k];
        }
    }

    inline void inject_fixed_params(
        const Eigen::VectorXd &x,
        std::vector<had::AReal> &p,
        const std::vector<int> &fixed_idx)
    {
        for (size_t k = 0; k < fixed_idx.size(); ++k)
            p[fixed_idx[k]] = had::AReal(x[k]);
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

    inline void inject_random_params(
        const std::vector<had::AReal> &u_ad,
        std::vector<had::AReal> &p,
        const std::vector<int> &random_idx)
    {
        for (size_t k = 0; k < random_idx.size(); ++k)
        {
            p[random_idx[k]] = u_ad[k];
        }
    }
    inline void inject_random_params(
        const std::vector<double> &u,
        std::vector<had::AReal> &p,
        const std::vector<int> &random_idx)
    {
        for (size_t k = 0; k < random_idx.size(); ++k)
        {
            p[random_idx[k]] = had::AReal(u[k]);
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

    using HessianPattern = std::vector<std::pair<int, int>>;

    inline std::unordered_map<int, HessianPattern> &pattern_cache()
    {
        static std::unordered_map<int, HessianPattern> cache;
        return cache;
    }

    inline HessianPattern discover_pattern_dense(
        const ADScope &scope,
        const std::vector<had::AReal> &p_full,
        const std::vector<int> &random_idx,
        double tol = 1e-12)
    {
        std::cout<<"Quadra: Discovering Hessian pattern (dense) for n_random = " << random_idx.size() << " ...\n";
        const int n = (int)random_idx.size();
        HessianPattern pattern;
        pattern.reserve(n * n);

        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                double hij = scope.hess(
                    p_full[random_idx[i]],
                    p_full[random_idx[j]]);
                if (std::abs(hij) > tol)
                    pattern.emplace_back(i, j);
            }
        }
        return pattern;
    }

    inline const HessianPattern &get_pattern(
        const ADScope &scope,
        const std::vector<had::AReal> &p_full,
        const std::vector<int> &random_idx)
    {
        const int n = (int)random_idx.size();
        auto &cache = pattern_cache();

        auto it = cache.find(n);
        if (it != cache.end())
            return it->second;

        auto pattern = discover_pattern_dense(scope, p_full, random_idx);
        auto res = cache.emplace(n, std::move(pattern));
        return res.first->second;
    }

    inline HessianPattern banded_hessian_pattern(int n, int bandwidth)
    {
        HessianPattern pattern;

        for (int i = 0; i < n; ++i)
        {
            int j0 = std::max(0, i - bandwidth);
            int j1 = std::min(n - 1, i + bandwidth);

            for (int j = j0; j <= j1; ++j)
                pattern.emplace_back(i, j);
        }

        return pattern;
    }

    inline HessianPattern dense_hessian_pattern(int n)
    {
        HessianPattern pattern;
        pattern.reserve(n * n);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                pattern.emplace_back(i, j);

        return pattern;
    }

    inline Eigen::SparseMatrix<double> extract_sparse_hessian(
        const ADScope &scope,
        const std::vector<had::AReal> &p_full,
        const std::vector<int> &random_idx,
        const HessianPattern &pattern,
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
            std::vector<had::AReal> p_full;
            p_full.reserve(params.size());

            for (int i = 0; i < params.size(); ++i)
            {
                p_full.emplace_back(had::AReal(0.0));
            }

            inject_fixed_params(x, p_full, fixed_idx);
            inject_random_params(u, p_full, random_idx);

            // --------------------------------------------------
            // Forward pass
            // --------------------------------------------------
            had::AReal nll = model(p_full);

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

        std::vector<had::AReal> p_full;
        p_full.reserve(params.size());

        for (int i = 0; i < params.size(); ++i)
        {
            p_full.emplace_back(had::AReal(0.0));
        }

        inject_fixed_params(x, p_full, fixed_idx);
        inject_random_params(u_star, p_full, random_idx);

        had::AReal nll = model(p_full);

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

        res.value = nll.val + 0.5 * logdet;

        return res;
    }

} // namespace quadra

#endif // LAPLACE_HPP