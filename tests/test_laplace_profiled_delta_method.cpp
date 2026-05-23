#include "../core/laplace/laplace_profiled_delta_method.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

int main()
{
    const double estimate = 10.0;

    Eigen::VectorXd g_theta(2);
    g_theta << 1.0, -2.0;

    Eigen::VectorXd g_u(2);
    g_u << 3.0, 4.0;

    Eigen::MatrixXd du_dtheta(2, 2);
    du_dtheta << 0.5, -0.1,
                 0.2,  0.3;

    Eigen::MatrixXd theta_covariance(2, 2);
    theta_covariance << 0.04, 0.01,
                        0.01, 0.09;

    const auto result =
        quadra::compute_laplace_profiled_delta_method(
            estimate,
            g_theta,
            g_u,
            du_dtheta,
            theta_covariance);

    if (!result.success_m) {
        std::cerr << "FAIL: profiled delta method failed: "
                  << result.message_m << "\n";
        return 1;
    }

    const Eigen::VectorXd expected_gradient =
        g_theta + (g_u.transpose() * du_dtheta).transpose();

    const double expected_variance =
        (expected_gradient.transpose() *
         theta_covariance *
         expected_gradient)(0, 0);

    const double gradient_error =
        (result.gradient_m - expected_gradient).norm();

    const double variance_error =
        std::abs(result.variance_m - expected_variance);

    if (gradient_error > 1.0e-12) {
        std::cerr << "FAIL: gradient mismatch\n";
        return 1;
    }

    if (variance_error > 1.0e-12) {
        std::cerr << "FAIL: variance mismatch\n";
        return 1;
    }

    if (std::abs(result.std_error_m - std::sqrt(expected_variance)) > 1.0e-12) {
        std::cerr << "FAIL: std.error mismatch\n";
        return 1;
    }

    if (std::abs(result.cv_m - result.std_error_m / estimate) > 1.0e-12) {
        std::cerr << "FAIL: CV mismatch\n";
        return 1;
    }

    std::cout << "PASS: Laplace profiled delta-method utility\n";
    std::cout << "  gradient: " << result.gradient_m.transpose() << "\n";
    std::cout << "  variance: " << result.variance_m << "\n";
    std::cout << "  std.error: " << result.std_error_m << "\n";
    std::cout << "  cv: " << result.cv_m << "\n";

    return 0;
}
