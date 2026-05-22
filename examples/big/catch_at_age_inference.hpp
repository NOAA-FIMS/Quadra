#ifndef QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP
#define QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP

#include "../../core/inference/fixed_effect_covariance.hpp"
#include "../../core/inference/fixed_effect_report.hpp"
#include "../../core/inference/delta_method.hpp"
#include "../../core/inference/ad_delta_method.hpp"
#include "../../core/inference/ad_delta_method_vector.hpp"
#include "../../core/laplace/laplace_profiled_derived_gradient.hpp"
#include "../../core/laplace/laplace_profiled_ad_gradient.hpp"
#include "catch_at_age_derived.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <type_traits>

namespace example {

template <typename Objective>
inline void run_big_catch_at_age_inference(
    Objective& objective,
    const CatchAtAgeLaplaceModel& model,
    const std::vector<double>& final_random_effects,
    const std::vector<std::string>& parameter_names,
    const std::vector<double>& theta_hat,
    const Eigen::MatrixXd& du_dtheta,
    bool du_dtheta_available)
{
    using namespace quadra;

    std::cout << "\n========================================\n";
    std::cout << "Quadra inference layer\n";
    std::cout << "========================================\n";

    auto covariance =
        estimate_fixed_effect_covariance(
            objective,
            theta_hat,
            1.0e-4);

    if (!covariance.success_m) {
        std::cout << "Covariance estimation failed: "
                  << covariance.message_m << "\n";
        return;
    }

    std::cout << "\nCovariance estimation succeeded\n";

    auto report =
        build_fixed_effect_report(
            parameter_names,
            theta_hat,
            covariance);

    print_fixed_effect_report(report);

    std::cout << "\nDerived quantity uncertainty\n";

    auto depletion_fn =
        [&model, &final_random_effects](const auto& theta) {
            return evaluate_terminal_depletion_ad(
                model,
                theta,
                final_random_effects);
        };

    auto ssb_fn =
        [&model, &final_random_effects](const auto& theta) {
            return evaluate_terminal_ssb_proxy_ad(
                model,
                theta,
                final_random_effects);
        };

    auto mean_f_fn =
        [&model, &final_random_effects](const auto& theta) {
            return evaluate_mean_f_ad(
                model,
                theta,
                final_random_effects);
        };

    auto depletion_dm =
        ad_delta_method_scalar(
            depletion_fn,
            theta_hat,
            covariance.covariance_m);

    auto ssb_dm =
        ad_delta_method_scalar(
            ssb_fn,
            theta_hat,
            covariance.covariance_m);

    auto mean_f_dm =
        ad_delta_method_scalar(
            mean_f_fn,
            theta_hat,
            covariance.covariance_m);

    auto print_dm =
        [](const std::string& name,
           const auto& x)
        {
            std::cout
                << "\n"
                << name
                << "\n";

            if (!x.success_m) {
                std::cout
                    << "  failed: "
                    << x.message_m
                    << "\n";
                return;
            }

            std::cout
                << "  estimate : " << x.estimate_m << "\n"
                << "  std.error: " << x.std_error_m << "\n"
                << "  lower95  : " << x.lower95_m << "\n"
                << "  upper95  : " << x.upper95_m << "\n";
        };

    print_dm("terminal depletion", depletion_dm);
    print_dm("terminal ssb proxy", ssb_dm);
    print_dm("mean fishing mortality", mean_f_dm);


    std::cout << std::endl;
    std::cout << "Joint derived quantity inference" << std::endl;

    auto joint_fn =
        [&model, &final_random_effects](const auto& theta)
    {
        using Scalar = std::decay_t<decltype(theta[0])>;
        return std::vector<Scalar>{
            evaluate_terminal_depletion_ad(
                model,
                theta,
                final_random_effects),

            evaluate_terminal_ssb_proxy_ad(
                model,
                theta,
                final_random_effects),

            evaluate_mean_f_ad(
                model,
                theta,
                final_random_effects)
        };
    };

    const auto joint_dm =
        quadra::ad_delta_method_vector(
            joint_fn,
            theta_hat,
            covariance.covariance_m);

    if (!joint_dm.success_m)
    {
        std::cout << "joint derived inference failed: "
                  << joint_dm.message_m << std::endl;
    }
    else
    {
        const std::vector<std::string> names{
            "depletion",
            "terminal_ssb",
            "mean_F"
        };

        std::cout << std::endl;
        std::cout << "Joint estimates and standard errors" << std::endl;

        std::cout
            << std::setw(20) << "quantity"
            << std::setw(18) << "estimate"
            << std::setw(18) << "std.error"
            << std::endl;

        for (Eigen::Index i = 0;
             i < joint_dm.estimate_m.size();
             ++i)
        {
            std::cout
                << std::setw(20) << names[static_cast<std::size_t>(i)]
                << std::setw(18) << joint_dm.estimate_m[i]
                << std::setw(18) << joint_dm.std_error_m[i]
                << std::endl;
        }

        std::cout << std::endl;
        std::cout << "Derived quantity covariance matrix" << std::endl;
        std::cout << joint_dm.covariance_m << std::endl;

        std::cout << std::endl;
        std::cout << "Derived quantity correlation matrix" << std::endl;
        std::cout << joint_dm.correlation_m << std::endl;

        std::cout << std::endl;
        std::cout << "Derived quantity Jacobian" << std::endl;
        std::cout << joint_dm.jacobian_m << std::endl;
    }

    if (du_dtheta_available)
    {
        std::cout << std::endl;
        std::cout << "Profiled derived quantity uncertainty" << std::endl;

        auto print_profiled =
            [&](const std::string& name,
                double estimate,
                const Eigen::VectorXd& g_theta,
                const Eigen::VectorXd& g_u)
        {
            const auto profiled =
                quadra::compute_laplace_profiled_derived_gradient(
                    g_theta,
                    g_u,
                    du_dtheta);

            if (!profiled.success_m)
            {
                std::cout << "  " << name
                          << ": FAILED: "
                          << profiled.message_m
                          << std::endl;
                return;
            }

            const double variance =
                quadra::delta_variance_from_gradient(
                    profiled.gradient_m,
                    covariance.covariance_m);

            const double se =
                variance > 0.0
                    ? std::sqrt(variance)
                    : std::numeric_limits<double>::quiet_NaN();

            const double cv =
                std::abs(estimate) > 1.0e-12
                    ? se / std::abs(estimate)
                    : std::numeric_limits<double>::quiet_NaN();

            std::cout
                << "  " << name
                << ": estimate=" << estimate
                << " se=" << se
                << " cv=" << cv
                << std::endl;
        };

        {
            auto grad_blocks =
                quadra::evaluate_profiled_ad_gradient_blocks(
                    [&model](const auto& theta, const auto& u) {
                        return evaluate_terminal_depletion_ad(
                            model,
                            theta,
                            std::vector<double>{});
                    },
                    theta_hat,
                    final_random_effects);

            if (!grad_blocks.success_m)
            {
                std::cout << "  terminal depletion: FAILED gradient blocks: "
                          << grad_blocks.message_m << std::endl;
            }
            else
            {
                print_profiled(
                    "terminal depletion",
                    grad_blocks.estimate_m,
                    grad_blocks.gradient_fixed_m,
                    grad_blocks.gradient_random_m);
            }
        }

        {
            auto grad_blocks =
                quadra::evaluate_profiled_ad_gradient_blocks(
                    [&model](const auto& theta, const auto& u) {
                        return evaluate_terminal_ssb_proxy_ad(
                            model,
                            theta,
                            std::vector<double>{});
                    },
                    theta_hat,
                    final_random_effects);

            if (!grad_blocks.success_m)
            {
                std::cout << "  terminal ssb proxy: FAILED gradient blocks: "
                          << grad_blocks.message_m << std::endl;
            }
            else
            {
                print_profiled(
                    "terminal ssb proxy",
                    grad_blocks.estimate_m,
                    grad_blocks.gradient_fixed_m,
                    grad_blocks.gradient_random_m);
            }
        }

        {
            auto grad_blocks =
                quadra::evaluate_profiled_ad_gradient_blocks(
                    [&model](const auto& theta, const auto& u) {
                        return evaluate_mean_f_ad(
                            model,
                            theta,
                            std::vector<double>{});
                    },
                    theta_hat,
                    final_random_effects);

            if (!grad_blocks.success_m)
            {
                std::cout << "  mean fishing mortality: FAILED gradient blocks: "
                          << grad_blocks.message_m << std::endl;
            }
            else
            {
                print_profiled(
                    "mean fishing mortality",
                    grad_blocks.estimate_m,
                    grad_blocks.gradient_fixed_m,
                    grad_blocks.gradient_random_m);
            }
        }
    }


    std::cout << std::endl;
    std::cout << "Laplace mode sensitivity availability" << std::endl;
    if (!du_dtheta_available)
    {
        std::cout << "  unavailable: implicit derivative solve failed" << std::endl;
    }
    else
    {
        std::cout << "  du_dtheta dims: "
                  << du_dtheta.rows()
                  << " x "
                  << du_dtheta.cols()
                  << std::endl;
        std::cout << "  du_dtheta max abs: "
                  << du_dtheta.cwiseAbs().maxCoeff()
                  << std::endl;
    }
}

} // namespace example

#endif
