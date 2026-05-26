#include "../../include/quadra/quadra.hpp"
#include "../../core/laplace/laplace_exact_objective_gradient.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH()

struct ReuseStateSpaceModel
{
    int n_state = 100;

    template <typename Context>
    void initialize(Context&) {}

    template <typename T, typename Context>
    T evaluate(const std::vector<T>& p, Context&) const
    {
        const T mu = p[0];
        const T log_sigma_obs = p[1];
        const T log_sigma_rw = p[2];

        const T sigma_obs = exp(log_sigma_obs);
        const T sigma_rw = exp(log_sigma_rw);

        T nll = T(0.0);

        for (int t = 0; t < n_state; ++t)
        {
            const T x_t = p[3 + t];

            const double y_t =
                1.0
                + 0.10 * std::sin(0.03 * static_cast<double>(t))
                + 0.05 * std::cos(0.11 * static_cast<double>(t));

            const T obs_resid =
                (T(y_t) - (mu + x_t)) / sigma_obs;

            nll += T(0.5) * obs_resid * obs_resid + log_sigma_obs;

            if (t == 0)
            {
                const T z = x_t / sigma_rw;
                nll += T(0.5) * z * z + log_sigma_rw;
            }
            else
            {
                const T x_prev = p[3 + t - 1];
                const T z = (x_t - x_prev) / sigma_rw;
                nll += T(0.5) * z * z + log_sigma_rw;
            }
        }

        return nll;
    }
};

template <typename Function>
double elapsed_ms(Function&& f)
{
    const auto start = std::chrono::steady_clock::now();
    f();
    const auto end = std::chrono::steady_clock::now();

    return std::chrono::duration<double, std::milli>(
        end - start).count();
}

int main()
{
    std::ofstream csv(
        "benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark.csv");

    if (!csv)
    {
        std::cerr << "failed to open reuse benchmark CSV\n";
        return 1;
    }

    csv
        << "iteration,total_ms,objective_ms,"
        << "tape_setup_ms,reverse_pass_ms,"
        << "gradient_extract_ms,total_gradient_ms,"
        << "laplace_objective,gradient_norm,success\n";

    ReuseStateSpaceModel model;

    quadra::ParameterSet parameters;

    parameters.add(
        "mu",
        1.0,
        quadra::ParameterTransform::Identity,
        false);

    parameters.add(
        "log_sigma_obs",
        std::log(0.4),
        quadra::ParameterTransform::Identity,
        false);

    parameters.add(
        "log_sigma_rw",
        std::log(0.3),
        quadra::ParameterTransform::Identity,
        false);

    std::vector<double> theta{
        1.0,
        std::log(0.4),
        std::log(0.3)
    };

    std::vector<double> random_initial;

    for (int t = 0; t < model.n_state; ++t)
    {
        parameters.add(
            "x_" + std::to_string(t),
            0.0,
            quadra::ParameterTransform::Identity,
            true);

        random_initial.push_back(0.0);
    }

    quadra::LaplaceExactObjectiveGradientOptions options;

    for (int iter = 0; iter < 20; ++iter)
    {
        theta[0] += 0.0001;

        quadra::LaplaceExactObjectiveGradientResult result;

        const double total_ms =
            elapsed_ms([&]() {
                result =
                    quadra::evaluate_laplace_exact_objective_gradient(
                        model,
                        theta,
                        random_initial,
                        parameters,
                        options);
            });

        const bool success =
            result.converged_m &&
            result.logdet_ok_m &&
            std::isfinite(result.gradient_norm_m);

        csv
            << iter << ","
            << total_ms << ","
            << result.objective_ms_m << ","
            << result.tape_setup_ms_m << ","
            << result.reverse_pass_ms_m << ","
            << result.gradient_extract_ms_m << ","
            << result.total_gradient_ms_m << ","
            << result.laplace_objective_m << ","
            << result.gradient_norm_m << ","
            << (success ? 1 : 0)
            << "\n";

        std::cout
            << "iter=" << iter
            << " total_ms=" << total_ms
            << " objective_ms=" << result.objective_ms_m
            << " reverse_ms=" << result.reverse_pass_ms_m
            << " success=" << success
            << "\n";

        if (!success)
        {
            std::cerr << "reuse benchmark failed\n";
            return 1;
        }
    }

    return 0;
}
