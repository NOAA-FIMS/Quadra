#include "../core/uncertainty/selected_inverse_diagonal.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  using Triplet = Eigen::Triplet<double>;

  // SPD tridiagonal matrix:
  // [ 2 -1  0 ]
  // [-1  2 -1 ]
  // [ 0 -1  2 ]
  //
  // inverse =
  // [0.75 0.50 0.25]
  // [0.50 1.00 0.50]
  // [0.25 0.50 0.75]
  Eigen::SparseMatrix<double> h(3, 3);
  std::vector<Triplet> triplets;
  triplets.emplace_back(0, 0, 2.0);
  triplets.emplace_back(1, 1, 2.0);
  triplets.emplace_back(2, 2, 2.0);
  triplets.emplace_back(0, 1, -1.0);
  triplets.emplace_back(1, 0, -1.0);
  triplets.emplace_back(1, 2, -1.0);
  triplets.emplace_back(2, 1, -1.0);
  h.setFromTriplets(triplets.begin(), triplets.end());

  const std::vector<int> indices = {0, 2};
  const auto result =
      quadra::uncertainty::selected_inverse_submatrix_from_spd_hessian(
          h, indices);

  if (!result.success) {
    std::cerr << "selected inverse submatrix failed: "
              << result.message << "\n";
    return 1;
  }

  Eigen::MatrixXd expected(2, 2);
  expected << 0.75, 0.25,
              0.25, 0.75;

  const double max_abs_diff =
      (result.covariance - expected).cwiseAbs().maxCoeff();

  if (max_abs_diff > 1.0e-10) {
    std::cerr << "submatrix mismatch; max_abs_diff = " << max_abs_diff << "\n";
    std::cerr << "got:\n" << result.covariance << "\n";
    std::cerr << "expected:\n" << expected << "\n";
    return 1;
  }

  std::cout << "selected_inverse_submatrix_test_passed\n";
  std::cout << result.covariance << "\n";

  return 0;
}
