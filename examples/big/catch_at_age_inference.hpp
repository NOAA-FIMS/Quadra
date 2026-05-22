#ifndef QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP
#define QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP

#include "../../core/inference/fixed_effect_covariance.hpp"
#include "../../core/inference/fixed_effect_report.hpp"
#include "../../core/inference/delta_method.hpp"
#include "catch_at_age_derived.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace example {

template <typename Objective>
inline void run_big_catch_at_age_inference(
    Objective& objective,
    const CatchAtAgeLaplaceModel& model,
    const std::vector<double>& final_random_effects,
    const std::vector<std::string>& parameter_names,
    const std::vector<double>& theta_hat)
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
        [&model, &final_random_effects](const std::vector<double>& theta) {
            return evaluate_catch_at_age_derived_quantities(
                model,
                theta,
                final_random_effects).terminal_depletion_m;
        };

    auto ssb_fn =
        [&model, &final_random_effects](const std::vector<double>& theta) {
            return evaluate_catch_at_age_derived_quantities(
                model,
                theta,
                final_random_effects).terminal_ssb_proxy_m;
        };

    auto mean_f_fn =
        [&model, &final_random_effects](const std::vector<double>& theta) {
            return evaluate_catch_at_age_derived_quantities(
                model,
                theta,
                final_random_effects).mean_f_m;
        };

    auto depletion_dm =
        delta_method_scalar(
            depletion_fn,
            theta_hat,
            covariance.covariance_m);

    auto ssb_dm =
        delta_method_scalar(
            ssb_fn,
            theta_hat,
            covariance.covariance_m);

    auto mean_f_dm =
        delta_method_scalar(
            mean_f_fn,
            theta_hat,
            covariance.covariance_m);

    auto print_dm =
        [](const std::string& name,
           const quadra::DeltaMethodResult& x)
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
}

} // namespace example

#endif
