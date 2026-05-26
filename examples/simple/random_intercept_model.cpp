#include "../../core/laplace/laplace_implicit_workspace.hpp"
#include "../../core/laplace/laplace_profiled_derived_report.hpp"
#include "../../core/inference/fixed_effect_covariance.hpp"
#include "../../core/inference/report_serialization.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

DECLARE_ADGRAPH()

struct RandomInterceptModel
{
    std::vector<double> y{1.2, 0.8, 1.5, 1.1, 0.9};

    template <typename Context>
    void initialize(Context&)
    {
    }

    template <typename T, typename Context>
    T evaluate(const std::vector<T>& p, Context&) const
    {
        return evaluate_impl<T>(p);
    }

    template <typename T>
    T evaluate(const std::vector<T>& p) const
    {
        return evaluate_impl<T>(p);
    }

    template <typename T>
    T evaluate_impl(const std::vector<T>& p) const
    {
        const T mu = p[0];
        const T log_sigma = p[1];
        const T u = p[2];

        const T sigma = exp(log_sigma);

        T nll = T(0.0);

        // Random intercept prior: u ~ N(0, 1).
        nll += T(0.5) * u * u;

        // Observation model: y_i ~ N(mu + u, sigma).
        for (double yi : y)
        {
            const T resid = T(yi) - (mu + u);
            nll += T(0.5) * (resid / sigma) * (resid / sigma)
                   + log_sigma;
        }

        // Weak priors keep the tiny example well-conditioned.
        nll += T(0.5) * (mu / T(10.0)) * (mu / T(10.0));
        nll += T(0.5) * (log_sigma / T(2.0)) * (log_sigma / T(2.0));

        return nll;
    }
};

int main()
{
    RandomInterceptModel model;

    quadra::ParameterSet parameters;
    parameters.add("mu", 0.0, quadra::ParameterTransform::Identity, false);
    parameters.add("log_sigma", std::log(0.5), quadra::ParameterTransform::Identity, false);
    parameters.add("u", 0.0, quadra::ParameterTransform::Identity, true);

    std::vector<double> theta_hat{0.0, std::log(0.5)};
    std::vector<double> random_initial{0.0};

    auto objective_for_covariance =
        [&model, &random_initial](const std::vector<double>& theta)
    {
        const auto partition =
            quadra::partition_parameters([&]() {
                quadra::ParameterSet p;
                p.add("mu", theta[0], quadra::ParameterTransform::Identity, false);
                p.add("log_sigma", theta[1], quadra::ParameterTransform::Identity, false);
                p.add("u", 0.0, quadra::ParameterTransform::Identity, true);
                return p;
            }());

        const auto laplace =
            quadra::evaluate_laplace_objective(
                model,
                theta,
                random_initial,
                partition);

        return laplace.laplace_objective_m;
    };

    const auto covariance =
        quadra::estimate_fixed_effect_covariance(
            objective_for_covariance,
            theta_hat);

    if (!covariance.success_m)
    {
        std::cerr << "fixed-effect covariance failed: "
                  << covariance.message_m << "\n";
        return 1;
    }

    const auto workspace =
        quadra::build_laplace_implicit_workspace(
            model,
            theta_hat,
            random_initial,
            parameters);

    if (!workspace.success_m)
    {
        std::cerr << "workspace failed: "
                  << workspace.message_m << "\n";
        return 1;
    }

    std::vector<quadra::ProfiledDerivedQuantity> quantities;

    quantities.push_back({
        "population_mean",
        [](const std::vector<quadra::AD>& theta,
           const std::vector<quadra::AD>& u)
        {
            return theta[0] + u[0];
        }
    });

    quantities.push_back({
        "sigma",
        [](const std::vector<quadra::AD>& theta,
           const std::vector<quadra::AD>&)
        {
            return exp(theta[1]);
        }
    });

    const auto report =
        quadra::compute_laplace_profiled_derived_report(
            quantities,
            theta_hat,
            workspace.u_hat_m,
            covariance.covariance_m,
            workspace);

    if (!report.success_m)
    {
        std::cerr << "profiled report failed: "
                  << report.message_m << "\n";
        return 1;
    }

    quadra::write_profiled_derived_report_csv(
        "examples/simple/random_intercept_report.csv",
        report);

    std::cout << "Random intercept profiled derived report\n";

    for (std::size_t i = 0; i < report.names_m.size(); ++i)
    {
        const Eigen::Index k = static_cast<Eigen::Index>(i);

        std::cout
            << "  " << report.names_m[i]
            << ": estimate=" << report.delta_m.estimate_m[k]
            << " se=" << report.delta_m.std_error_m[k]
            << " cv=" << report.delta_m.cv_m[k]
            << "\n";
    }

    std::cout << "  workspace total ms: "
              << workspace.total_ms_m
              << "\n";

    std::cout << "  wrote examples/simple/random_intercept_report.csv\n";

    return 0;
}
