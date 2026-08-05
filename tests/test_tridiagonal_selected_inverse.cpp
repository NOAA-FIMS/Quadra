#include "../core/laplace/sparse_huu_factorization.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

int main() {
  constexpr int n = 17;
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);
  Eigen::MatrixXd Hdot = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    H(i, i) = 2.5 + 0.03 * i;
    Hdot(i, i) = 0.2 - 0.01 * i;
    if (i + 1 < n) {
      H(i, i + 1) = H(i + 1, i) = -0.45;
      Hdot(i, i + 1) = Hdot(i + 1, i) = 0.07 + 0.002 * i;
    }
  }

  quadra::laplace::SparseHuuFactorization factorization(
      quadra::laplace::dense_to_sparse(H));
  const double actual =
      factorization.trace_inverse_times(quadra::laplace::dense_to_sparse(Hdot));
  const double expected = (H.inverse() * Hdot).trace();
  if (std::abs(actual - expected) > 1e-12) {
    std::cerr << "tridiagonal selected-inverse trace mismatch: got " << actual
              << ", expected " << expected << "\n";
    return 1;
  }

  std::cout << "PASS: tridiagonal selected-inverse trace\n";
  return 0;
}
