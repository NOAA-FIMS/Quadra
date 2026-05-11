#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP
#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <limits>
#include <string>
#include "../external/eigen/Eigen/Dense"
#include "../external/LBFGSpp/include/LBFGS.h"
#include "autodiff.hpp"
#include "laplace.hpp"

namespace quadra
{

    //==============================
    struct OptResult
    {
        std::vector<double> par;
        double value;
        int iterations;
    };

    inline Eigen::VectorXd to_eigen(const std::vector<double> &v)
    {
        Eigen::VectorXd x(v.size());
        for (size_t i = 0; i < v.size(); ++i)
            x(i) = v[i];
        return x;
    }

    //==============================
    // LBFGS Objective (VALUE ONLY)
    //==============================
    template <typename Model>
    class LBFGSObjective
    {
        int print_every = 1;
        inline void print(int iter, double fx, double gnorm)
        {
            std::cout << "L-BFGS: "
                      << "outer iter = " << std::setw(3) << iter
                      << ", fx = " << std::setw(14) << std::fixed << std::setprecision(6) << fx
                      << ", |grad| = " << std::setw(12) << std::fixed << std::setprecision(6) << gnorm
                      << "\n";
        }

    public:
        Model &model;
        ParameterVector &params;
        std::vector<int> fixed_idx;
        std::vector<int> random_idx;
        LaplaceOptions options;
        int iter = 0;
        double last_fx = std::numeric_limits<double>::quiet_NaN();
        Eigen::VectorXd last_grad;
        Eigen::VectorXd last_x;

        LBFGSObjective(Model &m,
                       ParameterVector &p,
                       std::vector<int> fixed,
                       std::vector<int> random,
                       const LaplaceOptions &opts = default_laplace_options())
            : model(m), params(p),
              fixed_idx(fixed), random_idx(random),
              options(opts)
        {
            laplace_pattern_cache().clear();
        }

        double operator()(const VectorXd &x, VectorXd &grad)
        {
            TapeContext tape;
            had::ADGraph &graph = tape.graph;

            iter++;

            // 1. solve random effects OUTSIDE Laplace
            std::vector<double> u_star = solve_random_effects_laplace(model, params, x,
                                                                      fixed_idx,
                                                                      random_idx, graph);

            // 2. single Laplace evaluation at fixed u*
            using Result = LaplaceResult<Model>;
            Result res;
            res = laplace_eval_at_u_star(
                model,
                params,
                fixed_idx,
                random_idx,
                x,
                u_star,
                graph,
                options);

            grad = to_eigen(res.grad_x);

            last_fx = res.value;
            last_grad = grad;
            last_x = x;

            if ((iter % print_every == 0) || iter == 1)
            {
                print(iter, res.value, grad.norm());
            }

            return res.value;
        }
    };

    //==============================
    // Optimize
    //==============================
    template <typename Model>
    OptResult optimize_lbfgs(
        Model &model,
        ParameterVector &params,
        const LaplaceOptions &options = default_laplace_options())
    {
        using namespace LBFGSpp;
        using namespace Eigen;

        auto fixed_idx = build_fixed_index(params);
        auto random_idx = build_random_index(params);

        if (fixed_idx.empty())
        {
            throw std::runtime_error(
                "No fixed parameters found — optimizer has zero dimension");
        }

        //--------------------------------------------------
        // Initial vector
        //--------------------------------------------------
        VectorXd x(fixed_idx.size());
        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            x[k] = params.params[fixed_idx[k]].value;
        }

        //--------------------------------------------------
        // Objective
        //--------------------------------------------------
        LBFGSObjective<Model> fun(model, params, fixed_idx, random_idx, options);

        //--------------------------------------------------
        // Solver
        //--------------------------------------------------
        LBFGSParam<double> param;
        param.max_iterations = 100;
        param.epsilon = 1e-6;

        LBFGSSolver<double> solver(param);

        double fx = std::numeric_limits<double>::quiet_NaN();
        int niter = 0;

        try
        {
            niter = solver.minimize(fun, x, fx);
        }
        catch (const std::runtime_error &e)
        {
            const double gnorm = fun.last_grad.norm();
            const double max_grad =
                (fun.last_grad.size() > 0)
                    ? fun.last_grad.cwiseAbs().maxCoeff()
                    : std::numeric_limits<double>::infinity();

            // LBFGS++ can throw a line-search failure after the objective is
            // already effectively converged, especially for large sparse
            // Laplace problems where function decrease is tiny relative to
            // the objective scale. Treat that case as a graceful convergence.
            const std::string msg = e.what();
            const bool line_search_failed =
                msg.find("line search") != std::string::npos ||
                msg.find("Line search") != std::string::npos;

            const double convergence_like_grad = 1e-3;

            if (line_search_failed && max_grad < convergence_like_grad)
            {
                std::cout
                    << "L-BFGS: line search failed after convergence-like gradient. (max|grad| = " << max_grad << ") "
                    << "\nL-BFGS: "
                    << "outer iter = " << std::setw(3) << niter
                    << ", fx = " << std::setw(14) << std::fixed << std::setprecision(6) << fx
                    << ", |grad| = " << std::setw(12) << std::fixed << std::setprecision(6) << gnorm
                    << std::endl;

                if (fun.last_x.size() == x.size())
                {
                    x = fun.last_x;
                }

                fx = fun.last_fx;
                niter = fun.iter;
            }
            else
            {
                throw;
            }
        }

        //--------------------------------------------------
        // Write back
        //--------------------------------------------------
        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            params.params[fixed_idx[k]].value = x[k];
        }

        return {
            std::vector<double>(x.data(), x.data() + x.size()),
            fx,
            niter};
    }

} // namespace quadra
#endif
