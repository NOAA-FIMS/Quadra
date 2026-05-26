#include "../core/laplace/laplace_mode_sensitivity.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

Eigen::Vector2d solve_mode(double theta) {
  Eigen::Matrix2d H;
  H << 1.5, 0.10, 0.10, 1.30;

  Eigen::Vector2d b;
  b << theta, -0.5 * theta;

  return H.ldlt().solve(b);
}

Eigen::Vector2d finite_difference_du_dtheta(double theta, double h = 1.0e-6) {
  return (solve_mode(theta + h) - solve_mode(theta - h)) / (2.0 * h);
}

int main() {
  const double theta = 2.0;

  Eigen::Matrix2d H_uu;
  H_uu << 1.5, 0.10, 0.10, 1.30;

  Eigen::Matrix<double, 2, 1> H_u_theta;
  H_u_theta << -1.0, 0.5;

  const auto result =
      quadra::solve_dense_laplace_mode_sensitivity(H_uu, H_u_theta);

  if (!result.success_m) {
    std::cerr << "FAIL: Laplace mode sensitivity failed: " << result.message_m
              << "\n";
    return 1;
  }

  const Eigen::Vector2d expected = finite_difference_du_dtheta(theta);

  const Eigen::Vector2d actual = result.du_dtheta_m.col(0);

  const double error = (actual - expected).norm();

  if (error > 1.0e-8) {
    std::cerr << "FAIL: Laplace mode sensitivity mismatch\n";
    std::cerr << "actual:\n" << actual << "\n";
    std::cerr << "expected:\n" << expected << "\n";
    std::cerr << "error: " << error << "\n";
    return 1;
  }

  std::cout << "PASS: dense Laplace mode sensitivity test\n";
  std::cout << "du/dtheta:\n" << actual << "\n";
  std::cout << "error: " << error << "\n";

  return 0;
}
