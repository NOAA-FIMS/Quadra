#include "../core/laplace/laplace_structure_report.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool ok, const std::string &message) {
  if (!ok) {
    throw std::runtime_error(message);
  }
}

void require_near(double x, double y, double tol, const std::string &message) {
  if (std::abs(x - y) > tol) {
    throw std::runtime_error(message + ": got " + std::to_string(x) +
                             ", expected " + std::to_string(y));
  }
}

void test_diagonal_report() {
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 3);
  H(0, 0) = 4.0;
  H(1, 1) = 9.0;
  H(2, 2) = 16.0;

  const auto report = quadra::summarize_laplace_hessian_structure(H, 1.0e-12);

  require(report.random_effects == 3, "diagonal: random_effects");
  require(report.total_entries == 9, "diagonal: total_entries");
  require(report.structural_nonzeros == 3, "diagonal: structural_nonzeros");
  require_near(report.structural_density, 3.0 / 9.0, 1.0e-12,
               "diagonal: structural_density");
  require(report.eigen_success, "diagonal: eigen_success");
  require(report.positive_definite, "diagonal: positive_definite");
  require_near(report.min_eigenvalue, 4.0, 1.0e-12, "diagonal: min_eigenvalue");
  require_near(report.max_eigenvalue, 16.0, 1.0e-12,
               "diagonal: max_eigenvalue");

  // Absolute curvature values are 16, 9, 4. 25 / 29 = 86.2%, so 90%
  // requires all three nonzero entries.
  bool found_90 = false;
  for (const auto &row : report.effective_sparsity) {
    if (row.label == "90%") {
      found_90 = true;
      require(row.entries_required == 3, "diagonal: 90% entries");
    }
  }
  require(found_90, "diagonal: found 90% row");

  // Diagonal-only matrix should have effective bandwidth zero for all
  // targets below 100%.
  for (const auto &row : report.effective_bandwidth) {
    if (row.label != "100%") {
      require(row.bandwidth == 0, "diagonal: non-100% bandwidth");
    }
  }
}

void test_tridiagonal_effective_bandwidth() {
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(4, 4);
  H.diagonal().array() = 10.0;

  for (int i = 0; i < 3; ++i) {
    H(i, i + 1) = -4.0;
    H(i + 1, i) = -4.0;
  }

  // Weak long-range tails. These make the matrix structurally denser than
  // tridiagonal, but most curvature remains in bandwidth 1.
  H(0, 2) = 0.1;
  H(2, 0) = 0.1;
  H(1, 3) = 0.1;
  H(3, 1) = 0.1;
  H(0, 3) = 0.01;
  H(3, 0) = 0.01;

  const auto report = quadra::summarize_laplace_hessian_structure(H, 1.0e-12);

  require(report.random_effects == 4, "tri: random_effects");
  require(report.total_entries == 16, "tri: total_entries");
  require(report.structural_nonzeros == 16, "tri: structural_nonzeros");
  require(report.positive_definite, "tri: positive_definite");

  std::size_t bw95 = 999;
  std::size_t entries95 = 0;
  for (const auto &row : report.effective_bandwidth) {
    if (row.label == "95%") {
      bw95 = row.bandwidth;
    }
  }
  for (const auto &row : report.effective_sparsity) {
    if (row.label == "95%") {
      entries95 = row.entries_required;
    }
  }

  require(bw95 == 1, "tri: 95% effective bandwidth should be 1");
  require(entries95 < report.structural_nonzeros,
          "tri: 95% effective sparsity should compress structural nnz");
}

void test_non_positive_definite_detection() {
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2, 2);
  H(0, 0) = 1.0;
  H(1, 1) = -1.0;

  const auto report = quadra::summarize_laplace_hessian_structure(H, 1.0e-12);

  require(report.eigen_success, "nonpd: eigen_success");
  require(!report.positive_definite, "nonpd: should not be positive definite");
  require(report.min_eigenvalue < 0.0, "nonpd: min eigenvalue negative");
}

} // namespace

int main() {
  try {
    test_diagonal_report();
    test_tridiagonal_effective_bandwidth();
    test_non_positive_definite_detection();
  } catch (const std::exception &e) {
    std::cerr << "test_laplace_structure_report failed: " << e.what() << "\n";
    return 1;
  }

  std::cout << "test_laplace_structure_report passed\n";
  return 0;
}
