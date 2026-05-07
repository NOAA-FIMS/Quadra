#ifndef QUADRA_COVARIANCE_HPP
#define QUADRA_COVARIANCE_HPP

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "../eigen/Eigen/Dense"
#include "../eigen/Eigen/Sparse"
#include "../eigen/Eigen/SparseCholesky"

#include "../../model/parameter.hpp"
#include "../laplace/laplace.hpp"
#include "../optimizer/optimizer.hpp"

namespace quadra
{

    struct CovarianceOptions
    {
        bool include_fixed_effects = true;
        bool include_random_effects = false;
        bool use_laplace_adjusted_hessian = true;

        // First-pass implementation:
        // central finite difference of optimized Laplace gradient.
        double fd_step = 1e-4;

        double jitter_initial = 1e-12;
        int jitter_max_attempts = 12;

        bool store_on_parameters = true;
        bool verbose = true;
    };

    struct CovarianceResult
    {
        Eigen::MatrixXd hessian;
        Eigen::MatrixXd covariance;
        std::vector<int> parameter_indices;
        bool success = false;
    };

    inline std::vector<int> covariance_fixed_indices(
        const ParameterVector &params)
    {
        std::vector<int> idx;

        for (int i = 0; i < params.size(); ++i)
        {
            const auto &p =
                params.params[static_cast<std::size_t>(i)];

            if (!p.is_random &&
                p.active &&
                p.estimate_covariance)
            {
                idx.push_back(i);
            }
        }

        return idx;
    }

    inline void store_standard_errors(
        ParameterVector &params,
        const CovarianceResult &result)
    {
        for (Eigen::Index i = 0; i < result.covariance.rows(); ++i)
        {
            const int pidx =
                result.parameter_indices[static_cast<std::size_t>(i)];

            auto &p =
                params.params[static_cast<std::size_t>(pidx)];

            p.covariance_index = static_cast<int>(i);

            const double v = result.covariance(i, i);

            p.std_error =
                (v >= 0.0 && std::isfinite(v))
                    ? std::sqrt(v)
                    : std::numeric_limits<double>::quiet_NaN();
        }
    }

    inline Eigen::MatrixXd invert_spd_with_adaptive_jitter(
        const Eigen::MatrixXd &H,
        const CovarianceOptions &options)
    {
        if (H.rows() != H.cols())
        {
            throw std::invalid_argument(
                "invert_spd_with_adaptive_jitter: H must be square");
        }

        Eigen::LLT<Eigen::MatrixXd> llt(H);

        if (llt.info() == Eigen::Success)
        {
            return llt.solve(
                Eigen::MatrixXd::Identity(H.rows(), H.cols()));
        }

        double jitter = options.jitter_initial;

        for (int attempt = 0; attempt < options.jitter_max_attempts; ++attempt)
        {
            Eigen::MatrixXd H_reg = H;
            H_reg.diagonal().array() += jitter;

            llt.compute(H_reg);

            if (llt.info() == Eigen::Success)
            {
                if (options.verbose)
                {
                    std::cout
                        << "Quadra covariance: Hessian inversion succeeded with jitter = "
                        << jitter
                        << "\n";
                }

                return llt.solve(
                    Eigen::MatrixXd::Identity(H.rows(), H.cols()));
            }

            jitter *= 10.0;
        }

        throw std::runtime_error(
            "invert_spd_with_adaptive_jitter: covariance Hessian not SPD");
    }

    template <typename Model>
    Eigen::VectorXd laplace_gradient_at_fixed(
        Model &model,
        ParameterVector &params,
        const Eigen::VectorXd &x,
        const LaplaceOptions &laplace_options)
    {
        const auto fixed_idx = build_fixed_index(params);
        const auto random_idx = build_random_index(params);

        had::ADGraph graph;

        std::vector<double> u_star =
            solve_random_effects_laplace(
                model,
                params,
                x,
                fixed_idx,
                random_idx,
                graph);

        auto res =
            laplace_eval_at_u_star(
                model,
                params,
                fixed_idx,
                random_idx,
                x,
                u_star,
                graph,
                laplace_options);

        return to_eigen(res.grad_x);
    }

    template <typename Model>
    Eigen::MatrixXd laplace_fixed_hessian_fd(
        Model &model,
        ParameterVector &params,
        const Eigen::VectorXd &x_hat,
        const LaplaceOptions &laplace_options,
        const CovarianceOptions &covariance_options)
    {
        const Eigen::Index n = x_hat.size();

        Eigen::MatrixXd H(n, n);
        H.setZero();

        const double eps = covariance_options.fd_step;

        for (Eigen::Index j = 0; j < n; ++j)
        {
            Eigen::VectorXd xp = x_hat;
            Eigen::VectorXd xm = x_hat;

            xp[j] += eps;
            xm[j] -= eps;

            Eigen::VectorXd gp =
                laplace_gradient_at_fixed(
                    model,
                    params,
                    xp,
                    laplace_options);

            Eigen::VectorXd gm =
                laplace_gradient_at_fixed(
                    model,
                    params,
                    xm,
                    laplace_options);

            H.col(j) =
                (gp - gm) / (2.0 * eps);
        }

        // Symmetrize to reduce finite-difference noise.
        H = 0.5 * (H + H.transpose());

        return H;
    }

    template <typename Model>
    CovarianceResult estimate_fixed_covariance(
        Model &model,
        ParameterVector &params,
        const std::vector<double> &theta_hat,
        const LaplaceOptions &laplace_options = default_laplace_options(),
        const CovarianceOptions &covariance_options = CovarianceOptions())
    {
        if (!covariance_options.include_fixed_effects)
        {
            throw std::runtime_error(
                "estimate_fixed_covariance: include_fixed_effects is false");
        }

        if (covariance_options.include_random_effects)
        {
            throw std::runtime_error(
                "estimate_fixed_covariance: random-effect covariance not implemented yet");
        }

        const auto fixed_idx_all =
            build_fixed_index(params);

        const auto cov_idx =
            covariance_fixed_indices(params);

        if (cov_idx.empty())
        {
            throw std::runtime_error(
                "estimate_fixed_covariance: no fixed parameters marked for covariance estimation");
        }

        if (theta_hat.size() != fixed_idx_all.size())
        {
            throw std::invalid_argument(
                "estimate_fixed_covariance: theta_hat size does not match number of active fixed effects");
        }

        Eigen::VectorXd x_hat(
            static_cast<Eigen::Index>(theta_hat.size()));

        for (std::size_t i = 0; i < theta_hat.size(); ++i)
        {
            x_hat[static_cast<Eigen::Index>(i)] = theta_hat[i];
        }

        // Synchronize optimized fixed values into params.
        for (std::size_t k = 0; k < fixed_idx_all.size(); ++k)
        {
            params.params[static_cast<std::size_t>(fixed_idx_all[k])].value =
                x_hat[static_cast<Eigen::Index>(k)];
        }

        Eigen::MatrixXd H_full =
            laplace_fixed_hessian_fd(
                model,
                params,
                x_hat,
                laplace_options,
                covariance_options);

        Eigen::MatrixXd H_cov(
            static_cast<Eigen::Index>(cov_idx.size()),
            static_cast<Eigen::Index>(cov_idx.size()));

        for (std::size_t i = 0; i < cov_idx.size(); ++i)
        {
            auto it_i =
                std::find(
                    fixed_idx_all.begin(),
                    fixed_idx_all.end(),
                    cov_idx[i]);

            if (it_i == fixed_idx_all.end())
            {
                throw std::runtime_error(
                    "estimate_fixed_covariance: covariance index not found in fixed index set");
            }

            const Eigen::Index ii =
                static_cast<Eigen::Index>(
                    std::distance(fixed_idx_all.begin(), it_i));

            for (std::size_t j = 0; j < cov_idx.size(); ++j)
            {
                auto it_j =
                    std::find(
                        fixed_idx_all.begin(),
                        fixed_idx_all.end(),
                        cov_idx[j]);

                if (it_j == fixed_idx_all.end())
                {
                    throw std::runtime_error(
                        "estimate_fixed_covariance: covariance index not found in fixed index set");
                }

                const Eigen::Index jj =
                    static_cast<Eigen::Index>(
                        std::distance(fixed_idx_all.begin(), it_j));

                H_cov(
                    static_cast<Eigen::Index>(i),
                    static_cast<Eigen::Index>(j)) =
                    H_full(ii, jj);
            }
        }

        Eigen::MatrixXd C =
            invert_spd_with_adaptive_jitter(
                H_cov,
                covariance_options);

        CovarianceResult result;
        result.hessian = H_cov;
        result.covariance = C;
        result.parameter_indices = cov_idx;
        result.success = true;

        if (covariance_options.store_on_parameters)
        {
            store_standard_errors(params, result);
        }

        if (covariance_options.verbose)
        {
            std::cout
                << "Quadra covariance: estimated fixed-effect covariance for "
                << cov_idx.size()
                << " parameters.\n";
        }

        return result;
    }

    template <typename Model>
    CovarianceResult estimate_fixed_covariance(
        Model &model,
        ParameterVector &params,
        const OptResult &fit,
        const LaplaceOptions &laplace_options = default_laplace_options(),
        const CovarianceOptions &covariance_options = CovarianceOptions())
    {
        return estimate_fixed_covariance(
            model,
            params,
            fit.par,
            laplace_options,
            covariance_options);
    }

} // namespace quadra

#endif
