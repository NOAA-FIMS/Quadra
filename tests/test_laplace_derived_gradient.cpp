#include "../core/laplace/laplace_derived_gradient.hpp"
#include "../core/laplace/laplace_mode_sensitivity.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

Eigen::Vector2d solve_mode(double theta)
{
    Eigen::Matrix2d H;

    H << 1.5, 0.10,
         0.10, 1.30;

    Eigen::Vector2d b;

    b << theta,
        -0.5 * theta;

    return H.ldlt().solve(b);
}

double derived_quantity(
    double theta,
    const Eigen::Vector2d& u)
{
    return
        theta * theta
        + 2.0 * u[0]
        - 3.0 * u[1];
}

double profiled_derived(double theta)
{
    return derived_quantity(
        theta,
        solve_mode(theta));
}

double finite_difference_profiled_gradient(
    double theta,
    double h = 1.0e-6)
{
    return
        (profiled_derived(theta + h)
         - profiled_derived(theta - h))
        / (2.0 * h);
}

int main()
{
    const double theta = 2.0;

    Eigen::Matrix2d H_uu;

    H_uu << 1.5, 0.10,
            0.10, 1.30;

    Eigen::Matrix<double, 2, 1> H_u_theta;

    H_u_theta << -1.0,
                  0.5;

    const auto mode_result =
        quadra::solve_dense_laplace_mode_sensitivity(
            H_uu,
            H_u_theta);

    if (!mode_result.success_m) {
        std::cerr << "FAIL: mode sensitivity failed\n";
        return 1;
    }

    Eigen::VectorXd g_theta(1);
    g_theta[0] = 2.0 * theta;

    Eigen::VectorXd g_u(2);
    g_u << 2.0,
           -3.0;

    const auto grad_result =
        quadra::compute_laplace_derived_gradient(
            g_theta,
            g_u,
            mode_result.du_dtheta_m);

    if (!grad_result.success_m) {
        std::cerr << "FAIL: derived gradient failed\n";
        return 1;
    }

    const double actual =
        grad_result.gradient_m[0];

    const double expected =
        finite_difference_profiled_gradient(theta);

    const double error =
        std::abs(actual - expected);

    if (!std::isfinite(error)) {
        std::cerr << "FAIL: non-finite derived gradient error\n";
        return 1;
    }

    if (error > 1.0e-8) {
        std::cerr << "FAIL: Laplace derived gradient mismatch\n";
        std::cerr << "actual: " << actual << "\n";
        std::cerr << "expected: " << expected << "\n";
        std::cerr << "error: " << error << "\n";
        return 1;
    }

    std::cout << "PASS: Laplace derived gradient validation\n";
    std::cout << "gradient: " << actual << "\n";
    std::cout << "FD gradient: " << expected << "\n";
    std::cout << "error: " << error << "\n";

    return 0;
}
