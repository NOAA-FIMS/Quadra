#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

template <typename Type>
Type random_intercept_joint(
    const std::vector<double>& y,
    const Type& mu,
    const Type& u
) {
    Type nll = Type(0.0);

    for (double yi : y) {
        Type r = Type(yi) - (mu + u);
        nll += Type(0.5) * r * r;
    }

    nll += Type(0.5) * u * u;

    return nll;
}

std::vector<double> simulate_data(size_t n) {
    std::mt19937 rng(1234);
    std::normal_distribution<double> dist(5.25, 1.0);

    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i) {
        y[i] = dist(rng);
    }

    return y;
}

template <typename F>
double time_ms(F&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

struct EvalResult {
    double f_m = 0.0;
    double grad_mu_m = 0.0;
    double grad_u_m = 0.0;
};

EvalResult rebuild_eval(
    const std::vector<double>& y,
    double mu_value,
    double u_value
) {
    quadra::TapeContext tape;
    quadra::ADScope scope(tape.graph);

    quadra::AD mu = mu_value;
    quadra::AD u = u_value;

    quadra::AD nll = random_intercept_joint(y, mu, u);

    scope.backward(nll);

    Eigen::VectorXd g =
        quadra::extract_gradient(std::vector<quadra::AD>{mu, u});

    return {quadra::value_of(nll), g[0], g[1]};
}

int main() {
    std::cout << "\nQuadra reusable tape random-intercept benchmark\n\n";

    std::cout
        << std::setw(8)  << "n"
        << std::setw(12) << "evals"
        << std::setw(18) << "rebuild ms"
        << std::setw(18) << "reuse ms"
        << std::setw(14) << "speedup"
        << std::setw(18) << "final grad_mu"
        << std::setw(18) << "final grad_u"
        << "\n";

    std::cout << std::string(106, '-') << "\n";

    const int n_evals = 100;

    for (size_t n : std::vector<size_t>{10, 100, 1000, 5000, 10000}) {
        std::vector<double> y = simulate_data(n);

        EvalResult rebuild_last;

        double rebuild_ms = time_ms([&]() {
            for (int i = 0; i < n_evals; ++i) {
                double t = static_cast<double>(i) / static_cast<double>(n_evals - 1);
                double mu = 4.0 + 1.0 * t;
                double u = 0.25 * (1.0 - t);

                rebuild_last = rebuild_eval(y, mu, u);
            }
        });

        EvalResult reuse_last;

        double reuse_ms = time_ms([&]() {
            quadra::TapeContext tape;
            quadra::ADScope scope(tape.graph);

            quadra::AD mu = 4.0;
            quadra::AD u = 0.25;

            quadra::AD nll = random_intercept_joint(y, mu, u);

            for (int i = 0; i < n_evals; ++i) {
                double t = static_cast<double>(i) / static_cast<double>(n_evals - 1);
                double mu_value = 4.0 + 1.0 * t;
                double u_value = 0.25 * (1.0 - t);

                quadra::set_value(mu, mu_value);
                quadra::set_value(u, u_value);

                scope.forward();
                scope.zero_adjoints();
                scope.backward(nll);

                Eigen::VectorXd g =
                    quadra::extract_gradient(std::vector<quadra::AD>{mu, u});

                reuse_last = {quadra::value_of(nll), g[0], g[1]};
            }
        });

        const double speedup = rebuild_ms / reuse_ms;

        std::cout
            << std::setw(8)  << n
            << std::setw(12) << n_evals
            << std::setw(18) << std::fixed << std::setprecision(3) << rebuild_ms
            << std::setw(18) << reuse_ms
            << std::setw(14) << speedup
            << std::setw(18) << reuse_last.grad_mu_m
            << std::setw(18) << reuse_last.grad_u_m
            << "\n";
    }

    return 0;
}
