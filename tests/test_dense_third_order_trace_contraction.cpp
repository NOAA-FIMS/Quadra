#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <vector>

#include "../core/laplace/sparse_huu_factorization.hpp"
#include "../core/laplace/third_order_dense_hdot.hpp"

DECLARE_ADGRAPH()

int main() {
  auto objective = [](const auto &x) {
    const auto scale = exp(x[0]);
    const auto sum = x[1] + 2.0 * x[2] - 0.5 * x[3];
    return 0.5 * (x[1] * x[1] + 2.0 * x[2] * x[2] + 3.0 * x[3] * x[3]) +
           0.5 * scale * sum * sum;
  };
  const std::vector<double> x = {0.2, 0.1, -0.2, 0.3};
  const std::vector<double> direction = {1.0, 0.03, -0.02, 0.01};
  const std::vector<int> random_indices = {1, 2, 3};

  Eigen::MatrixXd H(3, 3);
  const double s = std::exp(x[0]);
  Eigen::Vector3d loading(1.0, 2.0, -0.5);
  H.setZero();
  H.diagonal() << 1.0, 2.0, 3.0;
  H += s * loading * loading.transpose();
  quadra::laplace::SparseHuuFactorization factor(H.sparseView());

  const Eigen::MatrixXd hdot =
      quadra::laplace::dense_hdot_third_order_polarized(objective, x, direction,
                                                        random_indices);
  const double expected = (H.inverse() * hdot).trace();
  const double streamed =
      quadra::laplace::dense_hdot_trace_third_order_polarized(
          objective, x, direction, random_indices, [&](int column) {
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(3);
            rhs[column] = 1.0;
            return factor.solve(rhs);
          });

  if (std::abs(streamed - expected) > 1.0e-10) {
    std::cerr << "streamed dense trace mismatch: " << streamed << " vs "
              << expected << '\n';
    return 1;
  }
  Eigen::LLT<Eigen::MatrixXd> inverse_llt(H.inverse());
  const Eigen::MatrixXd inverse_factor = inverse_llt.matrixL();
  const double factorized =
      quadra::laplace::dense_hdot_trace_third_order_factorized(
          objective, x, direction, random_indices, inverse_factor);
  if (std::abs(factorized - expected) > 1.0e-10) {
    std::cerr << "factorized dense trace mismatch: " << factorized << " vs "
              << expected << '\n';
    return 1;
  }
  const double diagonal_only =
      quadra::laplace::dense_hdot_trace_third_order_polarized(
          objective, x, direction, random_indices,
          [&](int column) {
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(3);
            rhs[column] = 1.0;
            return factor.solve(rhs);
          },
          0);
  double expected_diagonal = 0.0;
  const Eigen::MatrixXd Hinv = H.inverse();
  for (int i = 0; i < 3; ++i)
    expected_diagonal += Hinv(i, i) * hdot(i, i);
  if (std::abs(diagonal_only - expected_diagonal) > 1.0e-10) {
    std::cerr << "band-limited dense trace mismatch\n";
    return 1;
  }
  std::cout << "streamed dense third-order trace test passed\n";
}
