#include "../core/laplace/sparse_factorization_cache.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  Eigen::SparseMatrix<double> H(3, 3);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.emplace_back(0, 0, 4.0);
  triplets.emplace_back(1, 1, 3.0);
  triplets.emplace_back(2, 2, 2.0);
  triplets.emplace_back(0, 1, 0.5);
  triplets.emplace_back(1, 0, 0.5);
  triplets.emplace_back(1, 2, 0.25);
  triplets.emplace_back(2, 1, 0.25);

  H.setFromTriplets(triplets.begin(), triplets.end());
  H.makeCompressed();

  Eigen::MatrixXd B(3, 2);
  B << 1.0, 0.5, 2.0, 1.0, -1.0, 3.0;

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> direct;
  direct.compute(H);

  if (direct.info() != Eigen::Success) {
    std::cerr << "FAIL: direct LDLT failed\n";
    return 1;
  }

  const Eigen::MatrixXd expected = direct.solve(B);

  quadra::SparseLDLTFactorizationCache cache;
  cache.analyze_pattern(H);
  cache.factorize(H);

  const Eigen::MatrixXd actual = cache.solve(B);

  const double expected_logdet = direct.vectorD().array().log().sum();
  const double logdet_error = std::abs(cache.logdet() - expected_logdet);

  const double error = (actual - expected).cwiseAbs().maxCoeff();

  if (!std::isfinite(error) || error > 1.0e-12) {
    std::cerr << "FAIL: cached LDLT solve mismatch: " << error << "\n";
    return 1;
  }

  if (!cache.analyzed() || !cache.factorized()) {
    std::cerr << "FAIL: cache flags not set\n";
    return 1;
  }

  cache.factorize(H);
  if (cache.symbolic_analysis_count() != 1 ||
      cache.numeric_factorization_count() != 2 || logdet_error > 1.0e-12) {
    std::cerr << "FAIL: sparse LDLT symbolic reuse/logdet mismatch\n";
    return 1;
  }

  std::cout << "PASS: sparse LDLT factorization cache\n";
  std::cout << "  error: " << error << "\n";
  std::cout << "  nnz: " << cache.nonzeros() << "\n";

  return 0;
}
