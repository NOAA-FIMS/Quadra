#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../core/laplace/full_exact_laplace_gradient_hdot.hpp"
#include "../core/laplace/had_quadra_hdot_provider.hpp"

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

void test_fd_hdot_provider_matches_analytic_hdot() {
  Eigen::VectorXd theta(3);
  theta << 1.1, std::log(0.45), std::log(0.8);

  Eigen::VectorXd uhat(5);
  uhat.setZero();

  auto hessian_uu = [](const Eigen::VectorXd &th, const Eigen::VectorXd &uh) {
    (void)uh;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(5, 5);
    const double h = std::exp(-2.0 * th[1]) + std::exp(-2.0 * th[2]);
    H.diagonal().array() = h;
    return H;
  };

  auto fd_hdot = quadra::laplace::make_finite_difference_hdot_provider(
      hessian_uu, 1.0e-6, true);

  for (int j = 0; j < 3; ++j) {
    const Eigen::MatrixXd got = fd_hdot(theta, uhat, j);

    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(5, 5);
    if (j == 1) {
      expected.diagonal().array() = -2.0 * std::exp(-2.0 * theta[1]);
    } else if (j == 2) {
      expected.diagonal().array() = -2.0 * std::exp(-2.0 * theta[2]);
    }

    const double max_abs = (got - expected).cwiseAbs().maxCoeff();
    std::cout << "theta_index " << j
              << " max |FD Hdot - analytic Hdot| = " << std::scientific
              << std::setprecision(10) << max_abs << "\n";

    if (max_abs > 1.0e-7) {
      throw std::runtime_error("finite-difference Hdot provider mismatch.");
    }
  }
}

void test_hdot_provider_plugs_into_gradient_layer() {
  Eigen::VectorXd theta(3);
  theta << 1.1, std::log(0.45), std::log(0.8);

  Eigen::VectorXd uhat(5);
  uhat.setZero();

  Eigen::VectorXd grad_joint(3);
  grad_joint << 0.25, 1.5, -0.7;

  auto hessian_uu = [](const Eigen::VectorXd &th, const Eigen::VectorXd &uh) {
    (void)uh;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(5, 5);
    const double h = std::exp(-2.0 * th[1]) + std::exp(-2.0 * th[2]);
    H.diagonal().array() = h;
    return H;
  };

  auto analytic_hdot = [](const Eigen::VectorXd &th, const Eigen::VectorXd &uh,
                          int j) {
    Eigen::MatrixXd Hdot = Eigen::MatrixXd::Zero(uh.size(), uh.size());
    if (j == 1) {
      Hdot.diagonal().array() = -2.0 * std::exp(-2.0 * th[1]);
    } else if (j == 2) {
      Hdot.diagonal().array() = -2.0 * std::exp(-2.0 * th[2]);
    }
    return Hdot;
  };

  auto fd_hdot = quadra::laplace::make_finite_difference_hdot_provider(
      hessian_uu, 1.0e-6, true);

  const Eigen::VectorXd grad_analytic =
      quadra::laplace::full_exact_laplace_gradient_hdot(
          grad_joint, hessian_uu, analytic_hdot, theta, uhat);

  const Eigen::VectorXd grad_fd_provider =
      quadra::laplace::full_exact_laplace_gradient_hdot(grad_joint, hessian_uu,
                                                        fd_hdot, theta, uhat);

  for (int j = 0; j < theta.size(); ++j) {
    expect_near(grad_fd_provider[j], grad_analytic[j], 1.0e-7,
                "FD Hdot provider plugged into gradient layer");
  }
}

} // namespace

int main() {
  test_fd_hdot_provider_matches_analytic_hdot();
  test_hdot_provider_plugs_into_gradient_layer();

  std::cout << "\nhad_quadra Hdot-provider scaffold tests passed\n";
  return 0;
}
