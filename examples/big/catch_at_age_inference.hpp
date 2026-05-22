#ifndef QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP
#define QUADRA_BIG_CATCH_AT_AGE_INFERENCE_HPP

#include "../../core/inference/fixed_effect_covariance.hpp"
#include "../../core/inference/fixed_effect_report.hpp"
#include "../../core/inference/delta_method.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace example {

template <typename Objective>
inline void run_big_catch_at_age_inference(
    Objective& objective,
    const std::vector<std::string>& parameter_names,
    const std::vector<double>& theta_hat,
    double terminal_ssb_proxy,
    double terminal_depletion,
    double mean_f)
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
        [terminal_depletion](const std::vector<double>&) {
            return terminal_depletion;
        };

    auto ssb_fn =
        [terminal_ssb_proxy](const std::vector<double>&) {
            return terminal_ssb_proxy;
        };

    auto mean_f_fn =
        [mean_f](const std::vector<double>&) {
            return mean_f;
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
