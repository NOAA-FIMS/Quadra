#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "../core/laplace/laplace_evaluator_exact_gradient_fd_adapter.hpp"

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

void test_adapter_matches_known_gradient() {
  Eigen::VectorXd theta(2);
  theta << 0.25, -0.6;

  Eigen::VectorXd uhat(2);
  uhat.setZero();

  auto joint_envelope_gradient = [](const Eigen::VectorXd &th,
                                    const Eigen::VectorXd & /*uhat*/) {
    Eigen::VectorXd g(2);
    g[0] = th[0];
    g[1] = 0.5 * th[1];
    return g;
  };

  auto hessian_uu_at_fixed_uhat = [](const Eigen::VectorXd &th,
                                     const Eigen::VectorXd & /*uhat*/) {
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 2);
    H(0, 0) = 1.0 + std::exp(th[0]);
    H(1, 1) = 2.0 + std::exp(th[1]);
    return H;
  };

  quadra::laplace::FullExactLaplaceGradientFDOptions options;
  options.step = 1.0e-6;
  options.relative_step = true;

  auto exact_gradient =
      quadra::laplace::make_laplace_evaluator_exact_gradient_fd_adapter(
          joint_envelope_gradient, hessian_uu_at_fixed_uhat, options);

  const Eigen::VectorXd got = exact_gradient(theta, uhat);

  Eigen::VectorXd expected(2);
  expected[0] =
      theta[0] + 0.5 * std::exp(theta[0]) / (1.0 + std::exp(theta[0]));
  expected[1] =
      0.5 * theta[1] + 0.5 * std::exp(theta[1]) / (2.0 + std::exp(theta[1]));

  expect_near(got[0], expected[0], 1.0e-8,
              "adapter full exact gradient theta[0]");
  expect_near(got[1], expected[1], 1.0e-8,
              "adapter full exact gradient theta[1]");
}

} // namespace

int main() {
  test_adapter_matches_known_gradient();
  std::cout << "laplace_evaluator_exact_gradient_fd_adapter tests passed\n";
  return 0;
}
