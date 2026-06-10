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
  // inverse diag = [0.75, 1.0, 0.75]
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

  const auto result =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(h);

  if (!result.success) {
    std::cerr << "selected inverse failed: " << result.message << "\n";
    return 1;
  }

  const double expected[] = {0.75, 1.0, 0.75};
  for (int i = 0; i < 3; ++i) {
    const double got = result.variance.at(static_cast<std::size_t>(i));
    if (std::abs(got - expected[i]) > 1.0e-10) {
      std::cerr << "variance mismatch at " << i << ": got " << got
                << ", expected " << expected[i] << "\n";
      return 1;
    }
  }

  std::cout << "selected_inverse_diagonal_test_passed\n";
  std::cout << "variance:";
  for (double v : result.variance)
    std::cout << " " << v;
  std::cout << "\nstandard_error:";
  for (double se : result.standard_error)
    std::cout << " " << se;
  std::cout << "\n";
  return 0;
}
