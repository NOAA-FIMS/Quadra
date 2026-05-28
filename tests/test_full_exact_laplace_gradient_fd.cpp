#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "../core/laplace/full_exact_laplace_gradient_fd.hpp"

namespace {

void expect_near(double got, double expected, double tol, const char *label) {
  const double err = std::abs(got - expected);
  if (!(err <= tol)) {
    std::cerr << "FAILED: " << label << "\n"
              << "  got      = " << got << "\n"
              << "  expected = " << expected << "\n"
              << "  abs err  = " << err << "\n"
              << "  tol      = " << tol << "\n";
    throw std::runtime_error(label);
  }
}

void test_scalar_theta_scalar_u() {
  // Joint objective:
  //
  //   f(theta, u) = 0.5 theta^2 + 0.5 (1 + exp(theta)) u^2
  //
  // Inner optimum:
  //
  //   uhat(theta) = 0
  //
  // Random-effect Hessian:
  //
  //   H_uu(theta, uhat) = 1 + exp(theta)
  //
  // Laplace objective up to constants:
  //
  //   L(theta) = 0.5 theta^2 + 0.5 log(1 + exp(theta))
  //
  // Full gradient:
  //
  //   dL/dtheta = theta + 0.5 exp(theta)/(1 + exp(theta))

  Eigen::VectorXd theta(1);
  theta[0] = 0.7;

  Eigen::VectorXd grad_joint_envelope(1);
  grad_joint_envelope[0] = theta[0];

  auto hessian_uu = [](const Eigen::VectorXd &th) {
    Eigen::MatrixXd H(1, 1);
    H(0, 0) = 1.0 + std::exp(th[0]);
    return H;
  };

  quadra::laplace::FullExactLaplaceGradientFDOptions options;
  options.step = 1.0e-6;
  options.relative_step = true;

  const Eigen::VectorXd got = quadra::laplace::full_exact_laplace_gradient_fd(
      grad_joint_envelope, hessian_uu, theta, options);

  const double expected =
      theta[0] + 0.5 * std::exp(theta[0]) / (1.0 + std::exp(theta[0]));

  expect_near(got[0], expected, 1.0e-8,
              "scalar theta / scalar random effect full Laplace gradient");
}

void test_two_theta_two_u_diagonal_hessian() {
  // H_uu(theta) = diag(1 + exp(theta_0), 2 + exp(theta_1))
  //
  // L(theta) = 0.5 theta.squaredNorm()
  //          + 0.5 log(1 + exp(theta_0))
  //          + 0.5 log(2 + exp(theta_1))
  //
  // This verifies the vector case and the trace contribution per parameter.

  Eigen::VectorXd theta(2);
  theta << -0.35, 0.9;

  Eigen::VectorXd grad_joint_envelope = theta;

  auto hessian_uu = [](const Eigen::VectorXd &th) {
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 2);
    H(0, 0) = 1.0 + std::exp(th[0]);
    H(1, 1) = 2.0 + std::exp(th[1]);
    return H;
  };

  quadra::laplace::FullExactLaplaceGradientFDOptions options;
  options.step = 1.0e-6;
  options.relative_step = true;

  const Eigen::VectorXd got = quadra::laplace::full_exact_laplace_gradient_fd(
      grad_joint_envelope, hessian_uu, theta, options);

  Eigen::VectorXd expected(2);
  expected[0] =
      theta[0] + 0.5 * std::exp(theta[0]) / (1.0 + std::exp(theta[0]));
  expected[1] =
      theta[1] + 0.5 * std::exp(theta[1]) / (2.0 + std::exp(theta[1]));

  expect_near(got[0], expected[0], 1.0e-8, "two theta / first trace component");
  expect_near(got[1], expected[1], 1.0e-8,
              "two theta / second trace component");
}

} // namespace

int main() {
  test_scalar_theta_scalar_u();
  test_two_theta_two_u_diagonal_hessian();

  std::cout << "full_exact_laplace_gradient_fd tests passed\n";
  return 0;
}
