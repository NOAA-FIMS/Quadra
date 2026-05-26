#include "../../include/quadra/quadra.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

DECLARE_ADGRAPH()

struct GaussianRandomWalkV1
{
    std::vector<double> y;

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
    T evaluate_impl(const std::vector<T>& p) const
    {
        const T mu = p[0];
        const T log_sigma_obs = p[1];
        const T log_sigma_rw = p[2];

        const T sigma_obs = exp(log_sigma_obs);
        const T sigma_rw = exp(log_sigma_rw);

        T nll = T(0.0);

        for (int t = 0; t < static_cast<int>(y.size()); ++t)
        {
            const T x_t = p[3 + t];

            const T obs_resid =
                (T(y[static_cast<std::size_t>(t)]) - (mu + x_t)) /
                sigma_obs;

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

        nll += T(0.001) * mu * mu;
        nll += T(0.001) * log_sigma_obs * log_sigma_obs;
        nll += T(0.001) * log_sigma_rw * log_sigma_rw;

        return nll;
    }
};

template <typename Function>
double elapsed_ms(Function&& f)
{
    const auto start = std::chrono::steady_clock::now();
    f();
    const auto end = std::chrono::steady_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

std::vector<double> make_observations(int n)
{
    std::vector<double> y;
    y.reserve(static_cast<std::size_t>(n));

    double x = 0.0;

    for (int t = 0; t < n; ++t)
    {
        x += 0.03 * std::sin(0.13 * static_cast<double>(t));
        const double obs =
            1.0 + x + 0.10 * std::cos(0.31 * static_cast<double>(t));

        y.push_back(obs);
    }

    return y;
}

int main()
{
    std::ofstream csv("benchmarks/state_space/state_space_benchmark.csv");

    if (!csv)
    {
        std::cerr << "failed to open benchmark CSV\n";
        return 1;
    }

    csv << "n_state,workspace_ms,implicit_derivatives_ms,factorization_ms,total_wall_ms,du_dtheta_max_abs,hessian_nnz,success\n";

    for (int n_state : {25, 50, 100, 250})
    {
        GaussianRandomWalkV1 model;
        model.y = make_observations(n_state);

        quadra::ParameterSet parameters;
        parameters.add("mu", 1.0, quadra::ParameterTransform::Identity, false);
        parameters.add("log_sigma_obs", std::log(0.4), quadra::ParameterTransform::Identity, false);
        parameters.add("log_sigma_rw", std::log(0.3), quadra::ParameterTransform::Identity, false);

        std::vector<double> theta{1.0, std::log(0.4), std::log(0.3)};
        std::vector<double> random_initial;
        random_initial.reserve(static_cast<std::size_t>(n_state));

        for (int t = 0; t < n_state; ++t)
        {
            parameters.add(
                "x_" + std::to_string(t),
                0.0,
                quadra::ParameterTransform::Identity,
                true);

            random_initial.push_back(0.0);
        }

        quadra::LaplaceImplicitWorkspace workspace;

        const double wall_ms =
            elapsed_ms([&]() {
                workspace =
                    quadra::build_laplace_implicit_workspace(
                        model,
                        theta,
                        random_initial,
                        parameters);
            });

        const double du_max =
            workspace.success_m
                ? workspace.du_dtheta_m.cwiseAbs().maxCoeff()
                : std::numeric_limits<double>::quiet_NaN();

        csv
            << n_state << ","
            << workspace.total_ms_m << ","
            << workspace.implicit_derivatives_ms_m << ","
            << workspace.factorization_ms_m << ","
            << wall_ms << ","
            << du_max << ","
            << workspace.H_uu_m.nonZeros() << ","
            << (workspace.success_m ? 1 : 0)
            << "\n";

        std::cout
            << "n_state=" << n_state
            << " workspace_ms=" << workspace.total_ms_m
            << " factorization_ms=" << workspace.factorization_ms_m
            << " nnz=" << workspace.H_uu_m.nonZeros()
            << " success=" << workspace.success_m
            << "\n";

        if (!workspace.success_m)
        {
            std::cerr << "workspace failed: "
                      << workspace.message_m
                      << "\n";
            return 1;
        }
    }

    std::cout << "wrote benchmarks/state_space/state_space_benchmark.csv\n";

    return 0;
}
