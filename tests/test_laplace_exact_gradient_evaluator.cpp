#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/laplace_exact_gradient_evaluator.hpp"

DECLARE_ADGRAPH()

namespace {

had::AReal objective(const std::vector<had::AReal>& theta,
                     const std::vector<had::AReal>& u) {
    const had::AReal& a = theta[0];
    const had::AReal& b = theta[1];

    had::AReal f(0.0);

    for (std::size_t i = 0; i < u.size(); ++i) {
        f = f + 0.5 * (u[i] - a) * (u[i] - a) + exp(b) * u[i] * u[i];
    }

    for (std::size_t i = 1; i < u.size(); ++i) {
        const had::AReal diff = u[i] - u[i - 1];
        f = f + 0.25 * diff * diff;
    }

    return f;
}

}  // namespace

int main() {
    using quadra::laplace::LaplaceExactGradientEvaluator;
    using quadra::laplace::MakeTridiagonalHdotPattern;

    std::vector<had::AReal> theta(2);
    std::vector<had::AReal> u(4);

    LaplaceExactGradientEvaluator evaluator;

    const Eigen::MatrixXd Hinv = Eigen::MatrixXd::Identity(4, 4);
    const auto pattern = MakeTridiagonalHdotPattern(4);
    const Eigen::VectorXd joint_gradient = Eigen::VectorXd::Zero(2);

    const auto result =
        evaluator.Evaluate(
            [&]() {
                theta[0] = had::AReal(0.2);
                theta[1] = had::AReal(-0.5);

                for (std::size_t i = 0; i < u.size(); ++i) {
                    u[i] = had::AReal(0.1 * static_cast<double>(i + 1));
                }

                return objective(theta, u);
            },
            &theta,
            &u,
            2,
            [](std::size_t k,
               Eigen::VectorXd& theta_direction,
               Eigen::VectorXd& random_direction) {
                theta_direction = Eigen::VectorXd::Zero(2);
                random_direction = Eigen::VectorXd::Zero(4);

                theta_direction[static_cast<int>(k)] = 1.0;

                for (int i = 0; i < random_direction.size(); ++i) {
                    random_direction[i] =
                        0.01 * static_cast<double>((i + 1) * (k + 1));
                }
            },
            [&](int row, int col) {
                return Hinv(row, col);
            },
            pattern,
            12.5,
            1.25,
            joint_gradient);

    if (std::abs(result.objective - 13.125) > 1.0e-12) {
        throw std::runtime_error("wrong Laplace exact-gradient objective");
    }

    if (result.gradient.size() != 2 || result.trace_terms.size() != 2) {
        throw std::runtime_error("wrong Laplace exact-gradient vector sizes");
    }

    if (!std::isfinite(result.gradient[0]) ||
        !std::isfinite(result.gradient[1])) {
        throw std::runtime_error("non-finite Laplace exact-gradient result");
    }

    std::cout << "laplace exact-gradient evaluator tests passed\n";
    std::cout << "objective = " << result.objective << "\n";
    std::cout << "gradient = " << result.gradient.transpose() << "\n";
    std::cout << "trace_terms = " << result.trace_terms.transpose() << "\n";

    return 0;
}
