#include <cmath>
#include <iostream>
#include <vector>

#include "../core/model/parameter.hpp"
#include "../core/autodiff/model_gradient.hpp"

DECLARE_ADGRAPH();

int main() {
    std::cout << "Testing transform log-Jacobian utilities\n";

    std::vector<double> x = {3.0, 0.25, 0.0, 2.0};

    std::vector<quadra::ParameterTransform> transforms = {
        quadra::ParameterTransform::Identity,
        quadra::ParameterTransform::Log,
        quadra::ParameterTransform::Logit,
        quadra::ParameterTransform::Square
    };

    auto terms = quadra::transform_log_jacobians(x, transforms);
    double total = quadra::sum_transform_log_jacobian(x, transforms);

    const double expected_identity = 0.0;
    const double expected_log = 0.25;
    const double expected_logit = std::log(0.5) + std::log(0.5);
    const double expected_square = std::log(std::abs(4.0));
    const double expected_total =
        expected_identity + expected_log + expected_logit + expected_square;

    std::cout << "identity logJ = " << terms[0] << "\n";
    std::cout << "log      logJ = " << terms[1] << "\n";
    std::cout << "logit    logJ = " << terms[2] << "\n";
    std::cout << "square   logJ = " << terms[3] << "\n";
    std::cout << "total    logJ = " << total << "\n";

    if (std::abs(terms[0] - expected_identity) > 1e-12) return 1;
    if (std::abs(terms[1] - expected_log) > 1e-12) return 1;
    if (std::abs(terms[2] - expected_logit) > 1e-12) return 1;
    if (std::abs(terms[3] - expected_square) > 1e-12) return 1;
    if (std::abs(total - expected_total) > 1e-12) return 1;

    // AD gradient check for log transform:
    // logJ = x, so d/dx logJ = 1.
    {
        quadra::TapeContext tape;
        quadra::ADScope scope(tape.graph);

        std::vector<quadra::AD> x_ad = quadra::to_ad(std::vector<double>{0.25});
        quadra::AD y =
            quadra::transform_log_jacobian(x_ad[0], quadra::ParameterTransform::Log);

        scope.backward(y);
        Eigen::VectorXd g = quadra::extract_gradient(x_ad);

        std::cout << "d logJ_log / dx = " << g[0] << "\n";

        if (std::abs(g[0] - 1.0) > 1e-12) {
            std::cerr << "FAIL: AD gradient for log transform Jacobian mismatch\n";
            return 1;
        }
    }

    // AD gradient check for logit transform at x = 0:
    // derivative is 1 - 2 * logistic(x), so zero at x = 0.
    {
        quadra::TapeContext tape;
        quadra::ADScope scope(tape.graph);

        std::vector<quadra::AD> x_ad = quadra::to_ad(std::vector<double>{0.0});
        quadra::AD y =
            quadra::transform_log_jacobian(x_ad[0], quadra::ParameterTransform::Logit);

        scope.backward(y);
        Eigen::VectorXd g = quadra::extract_gradient(x_ad);

        std::cout << "d logJ_logit / dx at 0 = " << g[0] << "\n";

        if (std::abs(g[0]) > 1e-12) {
            std::cerr << "FAIL: AD gradient for logit transform Jacobian mismatch\n";
            return 1;
        }
    }

    std::cout << "PASS\n";
    return 0;
}
