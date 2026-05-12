#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../core/laplace/laplace_exact_objective_lbfgs_optimizer.hpp"
#include "../core/laplace/laplace_objective.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class CorrelatedRandomInterceptModel :
    public quadra::QuadraModel<CorrelatedRandomInterceptModel> {
public:
    CorrelatedRandomInterceptModel(
        std::vector<double> y,
        std::vector<int> group,
        int n_groups,
        double rho,
        double lambda0,
        double lambda_diff
    )
        : y_m(std::move(y)),
          group_m(std::move(group)),
          n_groups_m(n_groups),
          rho_m(rho),
          lambda0_m(lambda0),
          lambda_diff_m(lambda_diff)
    {
        parameters_m.add(
            "mu",
            0.0,
            quadra::ParameterTransform::Identity,
            false
        );

        for (int g = 0; g < n_groups_m; ++g) {
            parameters_m.add(
                "u_" + std::to_string(g),
                0.0,
                quadra::ParameterTransform::Identity,
                true
            );
        }
    }

    std::vector<std::string> parameter_names_impl() const {
        return parameters_m.names();
    }

    const quadra::ParameterSet& parameters() const {
        return parameters_m;
    }

    template <typename Type>
    Type evaluate_impl(
        const std::vector<Type>& p,
        quadra::ModelReportContext&
    ) const {
        Type mu = p[0];

        Type nll = Type(0.0);

        for (size_t i = 0; i < y_m.size(); ++i) {
            const int g = group_m[i];
            Type u = p[1 + g];

            Type r = Type(y_m[i]) - (mu + u);
            nll += Type(0.5) * r * r;
        }

        Type u0 = p[1];
        nll += Type(0.5) * Type(lambda0_m) * u0 * u0;

        for (int g = 1; g < n_groups_m; ++g) {
            Type ug = p[1 + g];
            Type up = p[1 + g - 1];

            Type diff = ug - Type(rho_m) * up;
            nll += Type(0.5) * Type(lambda_diff_m) * diff * diff;
        }

        return nll;
    }

private:
    std::vector<double> y_m;
    std::vector<int> group_m;
    int n_groups_m;
    double rho_m;
    double lambda0_m;
    double lambda_diff_m;
    quadra::ParameterSet parameters_m;
};

struct SimData {
    std::vector<double> y_m;
    std::vector<int> group_m;
};

SimData simulate_grouped_data(
    int n_groups,
    int obs_per_group
) {
    std::mt19937 rng(1234);
    std::normal_distribution<double> u_dist(0.0, 0.25);
    std::normal_distribution<double> e_dist(0.0, 1.0);

    const double mu_true = 5.0;

    SimData out;
    out.y_m.reserve(static_cast<size_t>(n_groups * obs_per_group));
    out.group_m.reserve(static_cast<size_t>(n_groups * obs_per_group));

    double previous_u = 0.0;
    const double rho = 0.8;

    for (int g = 0; g < n_groups; ++g) {
        const double innovation = u_dist(rng);
        const double ug = rho * previous_u + innovation;
        previous_u = ug;

        for (int i = 0; i < obs_per_group; ++i) {
            out.y_m.push_back(mu_true + ug + e_dist(rng));
            out.group_m.push_back(g);
        }
    }

    return out;
}

template <typename F>
double time_ms(F&& f) {
    const auto start = std::chrono::high_resolution_clock::now();
    f();
    const auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

int main() {
    std::cout << "\nQuadra correlated-random-intercept scaling benchmark\n\n";

    std::cout
        << std::setw(8)  << "G"
        << std::setw(10) << "m"
        << std::setw(12) << "n"
        << std::setw(16) << "Laplace ms"
        << std::setw(16) << "LBFGS ms"
        << std::setw(12) << "iter"
        << std::setw(14) << "nnz Huu"
        << std::setw(16) << "theta_hat"
        << std::setw(16) << "grad_norm"
        << "\n";

    std::cout << std::string(120, '-') << "\n";

    const std::vector<int> groups = {5, 10, 25, 50, 100, 250};
    const int obs_per_group = 20;

    const double rho = 0.8;
    const double lambda0 = 1.0;
    const double lambda_diff = 1.0;

    for (int G : groups) {
        SimData data = simulate_grouped_data(G, obs_per_group);

        CorrelatedRandomInterceptModel model(
            data.y_m,
            data.group_m,
            G,
            rho,
            lambda0,
            lambda_diff
        );

        std::vector<double> theta0 = {4.0};
        std::vector<double> u0(static_cast<size_t>(G), 0.0);

        quadra::LaplaceObjectiveResult laplace_result;
        double laplace_ms = time_ms([&]() {
            laplace_result =
                quadra::evaluate_laplace_objective(
                    model,
                    theta0,
                    u0,
                    model.parameters()
                );
        });

        quadra::LaplaceExactObjectiveLBFGSResult opt_result;
        double opt_ms = time_ms([&]() {
            opt_result =
                quadra::optimize_laplace_fixed_effects_exact_objective_lbfgs(
                    model,
                    theta0,
                    u0,
                    model.parameters()
                );
        });

        const int n =
            static_cast<int>(data.y_m.size());

        std::cout
            << std::setw(8)  << G
            << std::setw(10) << obs_per_group
            << std::setw(12) << n
            << std::setw(16) << std::fixed << std::setprecision(3) << laplace_ms
            << std::setw(16) << opt_ms
            << std::setw(12) << opt_result.iterations_m
            << std::setw(14) << laplace_result.hessian_random_m.nonZeros()
            << std::setw(16) << opt_result.theta_hat_m[0]
            << std::setw(16) << opt_result.gradient_norm_m
            << "\n";
    }

    return 0;
}
