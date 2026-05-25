#ifndef QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP
#define QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP

#include "../../core/inference/fixed_effect_covariance.hpp"
#include "../../core/inference/fixed_effect_report.hpp"
#include "../../core/inference/delta_method.hpp"
#include "../../core/inference/ad_delta_method.hpp"
#include "../../core/inference/ad_delta_method_vector.hpp"
#include "../../core/laplace/laplace_profiled_derived_gradient.hpp"
#include "../../core/laplace/laplace_profiled_ad_gradient.hpp"
#include "../../core/laplace/laplace_profiled_delta_method.hpp"
#include "../../core/laplace/laplace_profiled_delta_method_vector.hpp"
#include "../../core/laplace/laplace_profiled_derived_report.hpp"
#include "../../core/laplace/laplace_implicit_workspace.hpp"
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
    const quadra::LaplaceImplicitWorkspace& implicit_workspace)
{
    using namespace quadra;

    std::cout << "\n========================================\n";
    std::cout << "Quadra inference layer\n";
    std::cout << "========================================\n";

    const bool du_dtheta_available =
        implicit_workspace.success_m;

    const Eigen::MatrixXd& du_dtheta =
        implicit_workspace.du_dtheta_m;

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
                quadra::compute_laplace_profiled_delta_method(
                    estimate,
                    g_theta,
                    g_u,
                    du_dtheta,
                    covariance.covariance_m);

            if (!profiled.success_m)
            {
                std::cout << "  " << name
                          << ": FAILED: "
                          << profiled.message_m
                          << std::endl;
                return;
            }

            const double fixed_u_variance =
                quadra::delta_variance_from_gradient(
                    g_theta,
                    covariance.covariance_m);

            const double fixed_u_se =
                fixed_u_variance > 0.0
                    ? std::sqrt(fixed_u_variance)
                    : std::numeric_limits<double>::quiet_NaN();

            const double se_abs_diff =
                profiled.std_error_m - fixed_u_se;

            const double se_rel_diff =
                std::abs(fixed_u_se) > 1.0e-12
                    ? se_abs_diff / fixed_u_se
                    : std::numeric_limits<double>::quiet_NaN();

            const double correction_norm =
                (profiled.gradient_m - g_theta).norm();

            std::cout
                << "  " << name
                << ": estimate=" << profiled.estimate_m
                << " fixed_u_se=" << fixed_u_se
                << " profiled_se=" << profiled.std_error_m
                << " profiled_cv=" << profiled.cv_m
                << " se_abs_diff=" << se_abs_diff
                << " se_rel_diff=" << se_rel_diff
                << " correction_norm=" << correction_norm
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



    if (du_dtheta_available)
    {
        std::cout << std::endl;
        std::cout << "Joint profiled derived quantity uncertainty" << std::endl;

        std::vector<quadra::ProfiledDerivedQuantity> quantities;

        quantities.push_back({
            "terminal_depletion",
            [&model](const std::vector<quadra::AD>& theta,
                     const std::vector<quadra::AD>& u)
            {
                return evaluate_terminal_depletion_ad(
                    model,
                    theta,
                    std::vector<double>{});
            }
        });

        quantities.push_back({
            "terminal_ssb_proxy",
            [&model](const std::vector<quadra::AD>& theta,
                     const std::vector<quadra::AD>& u)
            {
                return evaluate_terminal_ssb_proxy_ad(
                    model,
                    theta,
                    std::vector<double>{});
            }
        });

        quantities.push_back({
            "mean_fishing_mortality",
            [&model](const std::vector<quadra::AD>& theta,
                     const std::vector<quadra::AD>& u)
            {
                return evaluate_mean_f_ad(
                    model,
                    theta,
                    std::vector<double>{});
            }
        });

        const auto profiled_report =
            quadra::compute_laplace_profiled_derived_report(
                quantities,
                theta_hat,
                final_random_effects,
                covariance.covariance_m,
                implicit_workspace);

        if (!profiled_report.success_m)
        {
            std::cout << "  joint profiled derived uncertainty failed: "
                      << profiled_report.message_m
                      << std::endl;
        }
        else
        {
            std::cout
                << std::setw(28) << "quantity"
                << std::setw(18) << "estimate"
                << std::setw(18) << "std.error"
                << std::setw(18) << "cv"
                << std::endl;

            for (std::size_t i = 0;
                 i < profiled_report.names_m.size();
                 ++i)
            {
                std::cout
                    << std::setw(28) << profiled_report.names_m[i]
                    << std::setw(18) << profiled_report.delta_m.estimate_m[static_cast<Eigen::Index>(i)]
                    << std::setw(18) << profiled_report.delta_m.std_error_m[static_cast<Eigen::Index>(i)]
                    << std::setw(18) << profiled_report.delta_m.cv_m[static_cast<Eigen::Index>(i)]
                    << std::endl;
            }

            std::cout << std::endl;
            std::cout << "Profiled derived covariance matrix" << std::endl;
            std::cout << profiled_report.delta_m.covariance_m << std::endl;

            std::cout << std::endl;
            std::cout << "Profiled derived correlation matrix" << std::endl;
            std::cout << profiled_report.delta_m.correlation_m << std::endl;

            std::cout << std::endl;
            std::cout << "Profiled derived Jacobian" << std::endl;
            std::cout << profiled_report.delta_m.jacobian_m << std::endl;
        }
    }

    std::cout << std::endl;
    std::cout << "Laplace mode sensitivity availability" << std::endl;
    if (!du_dtheta_available)
    {
        std::cout << "  unavailable: "
                  << implicit_workspace.message_m
                  << std::endl;
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
