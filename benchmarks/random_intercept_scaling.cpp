#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_optimizer.hpp"
#include "../core/laplace/laplace_lbfgs_optimizer.hpp"
#include "../core/laplace/laplace_exact_gradient.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class RandomInterceptModel :
    public quadra::QuadraModel<RandomInterceptModel> {
public:
    RandomInterceptModel(std::vector<double> y)
        : y_m(std::move(y))
    {
        parameters_m.add(
            "mu",
            0.0,
            quadra::ParameterTransform::Identity,
            false
        );

        parameters_m.add(
            "u",
            0.0,
            quadra::ParameterTransform::Identity,
            true
        );
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
        Type u = p[1];

        Type nll = 0.0;

        for (double yi : y_m) {
            Type r = yi - (mu + u);
            nll += 0.5 * r * r;
        }

        nll += 0.5 * u * u;

        return nll;
    }

private:
    std::vector<double> y_m;
    quadra::ParameterSet parameters_m;
};

std::vector<double> simulate_data(
    size_t n,
    double mu,
    double u,
    double sigma,
    uint32_t seed
) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist(mu + u, sigma);

    std::vector<double> y(n);

    for (size_t i = 0; i < n; ++i) {
        y[i] = dist(rng);
    }

    return y;
}

double elapsed_ms(
    const std::chrono::high_resolution_clock::time_point& start,
    const std::chrono::high_resolution_clock::time_point& end
) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

int main() {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "Quadra random-intercept scaling benchmark\n";
    std::cout << "========================================\n\n";

    const std::vector<size_t> problem_sizes = {
        10,
        100,
        1000,
        5000,
        10000
    };

    std::cout
        << std::setw(8)  << "n"
        << std::setw(18) << "GD ms"
        << std::setw(14) << "GD iter"
        << std::setw(18) << "LBFGS ms"
        << std::setw(14) << "LBFGS iter"
        << std::setw(18) << "ExactGrad ms"
        << std::setw(18) << "theta_hat"
        << std::setw(18) << "u_hat"
        << "\n";

    std::cout << std::string(126, '-') << "\n";

    for (size_t n : problem_sizes) {
        std::vector<double> y =
            simulate_data(
                n,
                5.0,
                0.25,
                1.0,
                1234
            );

        RandomInterceptModel model(y);

        std::vector<double> theta0 = {4.0};
        std::vector<double> u0 = {0.0};

        quadra::LaplaceOptimizerOptions gd_options;
        gd_options.max_iterations_m = 200;
        gd_options.gradient_tolerance_m = 1e-6;
        gd_options.initial_step_scale_m = 0.5;

        quadra::LaplaceLBFGSOptions lbfgs_options;
        lbfgs_options.max_iterations_m = 200;
        lbfgs_options.memory_m = 5;
        lbfgs_options.gradient_tolerance_m = 1e-6;

        auto t0 = std::chrono::high_resolution_clock::now();

        auto gd_result =
            quadra::optimize_laplace_fixed_effects(
                model,
                theta0,
                u0,
                model.parameters(),
                gd_options
            );

        auto t1 = std::chrono::high_resolution_clock::now();

        auto lbfgs_result =
            quadra::optimize_laplace_fixed_effects_lbfgs(
                model,
                theta0,
                u0,
                model.parameters(),
                lbfgs_options
            );

        auto t2 = std::chrono::high_resolution_clock::now();

        auto exact_result =
            quadra::evaluate_laplace_exact_gradient(
                model,
                {lbfgs_result.theta_hat_m[0]},
                {lbfgs_result.u_hat_m[0]},
                model.parameters()
            );

        auto t3 = std::chrono::high_resolution_clock::now();

        const double gd_ms = elapsed_ms(t0, t1);
        const double lbfgs_ms = elapsed_ms(t1, t2);
        const double exact_ms = elapsed_ms(t2, t3);

        std::cout
            << std::setw(8)  << n
            << std::setw(18) << std::fixed << std::setprecision(3) << gd_ms
            << std::setw(14) << gd_result.iterations_m
            << std::setw(18) << lbfgs_ms
            << std::setw(14) << lbfgs_result.iterations_m
            << std::setw(18) << exact_ms
            << std::setw(18) << lbfgs_result.theta_hat_m[0]
            << std::setw(18) << lbfgs_result.u_hat_m[0]
            << "\n";
    }

    std::cout << "\n";
    std::cout << "Benchmark complete.\n";

    return 0;
}
