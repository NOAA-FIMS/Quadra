#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP
#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include "eigen/Eigen/Dense"
#include "LBFGSpp/include/LBFGS.h"
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
            std::cout << "[LBFGS] "
                      << "iter = " << std::setw(3) << iter
                      << ", fx = " << std::setw(14) << std::fixed << std::setprecision(6) << fx
                      << ", |grad| = " << std::setw(12) << std::fixed << std::setprecision(6) << gnorm
                      << "\n";
        }

    public:
        Model &model;
        ParameterVector &params;
        std::vector<int> fixed_idx;
        std::vector<int> random_idx;
        int iter = 0;

        LBFGSObjective(Model &m,
                       ParameterVector &p,
                       std::vector<int> fixed,
                       std::vector<int> random)
            : model(m), params(p),
              fixed_idx(fixed), random_idx(random)
        {
            pattern_cache().clear();
        }

        double operator()(const VectorXd &x, VectorXd &grad)
        {
            had::ADGraph graph;
            graph.Clear();
            had::g_ADGraph = &graph;

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
                graph);

            grad = to_eigen(res.grad_x);

            if ((iter % print_every == 0) || iter == 1)
            {
                print(iter, res.value, grad.norm());
            }

            return res.value;

            //    auto res = laplace_eval(
            //         model,
            //         params,
            //         random_idx,
            //         fixed_idx,
            //         x
            //     );

            //     grad = res.grad_x;
            //     return res.value;
        }
    };

    //==============================
    // Optimize
    //==============================
    template <typename Model>
    OptResult optimize_lbfgs(
        Model &model,
        ParameterVector &params)
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
        LBFGSObjective<Model> fun(model, params, fixed_idx, random_idx);

        //--------------------------------------------------
        // Solver
        //--------------------------------------------------
        LBFGSParam<double> param;
        param.max_iterations = 100;
        param.epsilon = 1e-6;

        LBFGSSolver<double> solver(param);

        double fx;
        int niter = solver.minimize(fun, x, fx);

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