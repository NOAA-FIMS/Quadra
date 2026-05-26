#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP
#pragma once

#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../external/eigen/Eigen/Dense"
#include "../external/LBFGSpp/include/LBFGS.h"

#include "autodiff.hpp"
#include "laplace.hpp"

namespace quadra
{

    struct OptResult
    {
        std::vector<double> par;
        double value = std::numeric_limits<double>::quiet_NaN();
        int iterations = 0;
        double grad_norm = std::numeric_limits<double>::quiet_NaN();
    };

    inline Eigen::VectorXd to_eigen(const std::vector<double> &x)
    {
        Eigen::VectorXd out(static_cast<Eigen::Index>(x.size()));
        for (Eigen::Index i = 0; i < out.size(); ++i)
        {
            out[i] = x[static_cast<size_t>(i)];
        }
        return out;
    }

    inline bool all_finite_eigen(const Eigen::VectorXd &v)
    {
        for (Eigen::Index i = 0; i < v.size(); ++i)
        {
            if (!std::isfinite(v[i]))
            {
                return false;
            }
        }
        return true;
    }

    inline double safe_eigen_norm(const Eigen::VectorXd &v)
    {
        if (!all_finite_eigen(v))
        {
            return std::numeric_limits<double>::infinity();
        }
        return v.norm();
    }

    template <typename Model>
    class LBFGSObjective
    {
        void print(int iter, double fx, double gnorm)
        {
            const bool converged = std::isfinite(gnorm) && gnorm <= epsilon;
            std::cout << "L-BFGS: "
                      << "outer eval = " << std::setw(3) << iter
                      << ", fx = " << std::setw(14) << std::fixed << std::setprecision(6) << fx
                      << ", |grad| = ";

            if (converged)
            {
                std::cout << "\033[1;32m";
            }
            else
            {
                std::cout << "\033[1;31m";
            }

            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << gnorm
                      << "\033[0m" << std::endl;
        }

    public:
        double epsilon = 1e-6;
        Model &model;
        ParameterVector &params;
        std::vector<int> fixed_idx;
        std::vector<int> random_idx;
        LaplaceOptions options;

        int iter = 0;
        int print_every = 10;

        double last_fx = std::numeric_limits<double>::quiet_NaN();
        Eigen::VectorXd last_grad;
        Eigen::VectorXd last_x;

        LBFGSObjective(
            Model &m,
            ParameterVector &p,
            std::vector<int> fixed,
            std::vector<int> random,
            const LaplaceOptions &opts = default_laplace_options())
            : model(m),
              params(p),
              fixed_idx(std::move(fixed)),
              random_idx(std::move(random)),
              options(opts)
        {
            laplace_pattern_cache().clear();
        }

        double operator()(const VectorXd &x, VectorXd &grad)
        {
            TapeContext tape;
            had::ADGraph &graph = tape.graph;

            ++iter;

            std::vector<double> u_star;
            const bool verbose_inner = ((iter % print_every) == 0) || iter == 1;

            try
            {
                u_star = solve_random_effects_laplace(
                    model,
                    params,
                    x,
                    fixed_idx,
                    random_idx,
                    graph);
            }
            catch (const std::exception &e)
            {
                std::cerr << "L-BFGS: random-effect mode solve failed; returning penalty. reason="
                          << e.what() << std::endl;
                const double penalty_gradient_scale = 1.0e3;

                for (int i = 0; i < grad.size(); ++i) {
                    const double xi = (i < x.size() && std::isfinite(x[i])) ? x[i] : 1.0;
                    grad[i] = penalty_gradient_scale * ((xi == 0.0) ? 1.0 : xi);
                }

                return std::numeric_limits<double>::max() / 100.0;
            }

            using Result = LaplaceResult<Model>;
            Result res;

            try
            {
                res = laplace_eval_at_u_star(
                    model,
                    params,
                    fixed_idx,
                    random_idx,
                    x,
                    u_star,
                    graph,
                    options);
            }
            catch (const std::exception &e)
            {
                std::cerr
                    << "L-BFGS: Laplace evaluation failed; returning penalty. reason="
                    << e.what()
                    << std::endl;

                grad.resize(x.size());
                grad.setZero();
                last_grad = grad;
                last_x = x;
                last_fx = std::numeric_limits<double>::max() / 1.0e100;
                return last_fx;
            }
            catch (...)
            {
                std::cerr
                    << "L-BFGS: Laplace evaluation failed with unknown exception; returning penalty."
                    << std::endl;

                grad.resize(x.size());
                grad.setZero();
                last_grad = grad;
                last_x = x;
                last_fx = std::numeric_limits<double>::max() / 1.0e100;
                return last_fx;
            }

            grad = to_eigen(res.grad_x);

            last_fx = res.value;
            last_grad = grad;
            last_x = x;

            const double gnorm = safe_eigen_norm(grad);

            if (gnorm <= epsilon || (iter % print_every) == 0 || iter == 1)
            {
                print(iter, res.value, gnorm);
            }

            return res.value;
        }
    };

    template <typename Model>
    OptResult optimize_lbfgs(
        Model &model,
        ParameterVector &params,
        const LaplaceOptions &options = default_laplace_options())
    {
        using namespace LBFGSpp;
        using namespace Eigen;

        const auto fixed_idx = build_fixed_index(params);
        const auto random_idx = build_random_index(params);

        if (fixed_idx.empty())
        {
            throw std::runtime_error(
                "No fixed parameters found — optimizer has zero dimension");
        }

        VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            x[static_cast<Eigen::Index>(k)] =
                params.params[static_cast<size_t>(fixed_idx[k])].value;
        }

        LBFGSObjective<Model> fun(model, params, fixed_idx, random_idx, options);
        fun.print_every = 10;

        LBFGSParam<double> param;
        param.max_iterations = 400;
        // param.max_linesearch = 20;
        param.epsilon = 1.0e-4;
        fun.epsilon = param.epsilon;

        LBFGSSolver<double> solver(param);

        double fx = std::numeric_limits<double>::quiet_NaN();
        int niter = 0;

        try
        {
            niter = solver.minimize(fun, x, fx);

            // quadra_lbfgs_honest_convergence_report_v1
            double quadra_final_fixed_grad_norm = std::numeric_limits<double>::quiet_NaN();
            if (fun.last_grad.size() > 0)
            {
                quadra_final_fixed_grad_norm = 0.0;
                for (int quadra_i = 0; quadra_i < fun.last_grad.size(); ++quadra_i)
                {
                    quadra_final_fixed_grad_norm += fun.last_grad[quadra_i] * fun.last_grad[quadra_i];
                }
                quadra_final_fixed_grad_norm = std::sqrt(quadra_final_fixed_grad_norm);
            }

            const bool quadra_requested_tol_met =
                std::isfinite(quadra_final_fixed_grad_norm) &&
                quadra_final_fixed_grad_norm <= 1.0e-4;

            std::cout << "L-BFGS minimize status report" << std::endl;
            std::cout << "  iterations returned by solver: " << niter << std::endl;
            std::cout << "  final objective returned by solver: " << fx << std::endl;
            std::cout << "  final fixed-gradient norm: " << quadra_final_fixed_grad_norm << std::endl;
            std::cout << "  requested gradient tolerance: " << std::scientific << 1.0e-4 << std::defaultfloat << std::endl;
            std::cout << "  configured max-iteration field: " << 400 << " (LBFGSpp max_iterations)" << std::endl;
            std::cout << "  requested tolerance met: "
                      << (quadra_requested_tol_met ? "yes" : "no") << std::endl;
            std::cout << "  outer convergence interpretation: "
                      << (quadra_requested_tol_met
                              ? "converged to requested gradient tolerance"
                              : "stopped before requested gradient tolerance; inspect LBFGS status/max iterations/line search")
                      << std::endl;
        }
        catch (const std::runtime_error &e)
        {
            const double gnorm = safe_eigen_norm(fun.last_grad);
            const double max_grad =
                (fun.last_grad.size() > 0)
                    ? fun.last_grad.cwiseAbs().maxCoeff()
                    : std::numeric_limits<double>::infinity();

            const std::string msg = e.what();
            const bool line_search_failed =
                msg.find("line search") != std::string::npos ||
                msg.find("Line search") != std::string::npos;

            const double convergence_like_grad = 2e-3;

            if (gnorm <= param.epsilon)
            {
                std::cout
                    << "L-BFGS: optimization reached convergence criterion "
                    << "(|grad| <= epsilon). max|grad| = " << max_grad
                    << std::endl;

                if (fun.last_x.size() == x.size())
                {
                    x = fun.last_x;
                }

                fx = fun.last_fx;
                niter = fun.iter;
            }
            else if (line_search_failed && max_grad < convergence_like_grad)
            {
                std::cout
                    << "L-BFGS: line search failed after convergence-like gradient. "
                    << "max|grad| = " << max_grad
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

        for (size_t k = 0; k < fixed_idx.size(); ++k)
        {
            params.params[static_cast<size_t>(fixed_idx[k])].value =
                x[static_cast<Eigen::Index>(k)];
        }

        OptResult result;

        if (fun.last_x.size() == x.size())
        {
            result.par.assign(fun.last_x.data(), fun.last_x.data() + fun.last_x.size());
        }
        else
        {
            result.par.assign(x.data(), x.data() + x.size());
        }

        result.value = std::isfinite(fun.last_fx) ? fun.last_fx : fx;
        result.iterations = niter;

        const double final_grad_norm = safe_eigen_norm(fun.last_grad);
        result.grad_norm =
            std::isfinite(final_grad_norm)
                ? final_grad_norm
                : std::numeric_limits<double>::infinity();

        return result;
    }

} // namespace quadra

#endif
