#include "../core/uncertainty/selected_inverse_diagonal.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <iostream>
#include <utility>
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

  std::vector<std::pair<int, int>> pairs = {
      {0, 0}, {1, 0}, {2, 0}, {2, 1}, {2, 2}};

  const auto result =
      quadra::uncertainty::selected_inverse_entries_from_spd_hessian(h, pairs);

  if (!result.success) {
    std::cerr << "selected inverse entries failed: " << result.message << "\n";
    return 1;
  }

  const double expected[] = {0.75, 0.50, 0.25, 0.50, 0.75};
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const double got = result.entries[i].covariance;
    if (std::abs(got - expected[i]) > 1.0e-10) {
      std::cerr << "entry mismatch at " << i << ": got " << got << ", expected "
                << expected[i] << "\n";
      return 1;
    }
  }

  std::cout << "selected_inverse_entries_test_passed\n";
  for (const auto &e : result.entries) {
    std::cout << e.row << "," << e.col << "," << e.covariance << "\n";
  }

  return 0;
}
