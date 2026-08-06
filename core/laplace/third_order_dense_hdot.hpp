#pragma once

#include "../../external/eigen/Eigen/Dense"
#include "../had_quadra.hpp"

#include <stdexcept>
#include <vector>

namespace quadra {
namespace laplace {

// Exact mixed third derivatives from the cubic directional derivative.
// The polarization identity extracts D^3 f[total_direction, e_a, e_b]
// without differencing Hessians or relying on reverse edge-pushing.
template <class Objective>
Eigen::MatrixXd
dense_hdot_third_order_polarized(Objective &&objective,
                                 const std::vector<double> &x,
                                 const std::vector<double> &total_direction,
                                 const std::vector<int> &random_indices) {
  if (x.size() != total_direction.size())
    throw std::invalid_argument("dense Hdot direction has wrong length");

  const int nr = static_cast<int>(random_indices.size());
  Eigen::MatrixXd out(nr, nr);
  for (int a = 0; a < nr; ++a) {
    for (int b = 0; b <= a; ++b) {
      auto cubic = [&](int st, int sa, int sb) {
        std::vector<double> direction = total_direction;
        for (double &value : direction)
          value *= st;
        direction[static_cast<std::size_t>(random_indices[a])] += sa;
        direction[static_cast<std::size_t>(random_indices[b])] += sb;
        return had::third_directional_derivative(objective, x, direction);
      };
      // Four-term polarization for a symmetric trilinear form. The usual
      // eight-term identity pairs opposite directions, whose cubic
      // directional derivatives differ only by sign.
      const double mixed =
          cubic(1, 1, 1) - cubic(1, 1, -1) - cubic(1, -1, 1) - cubic(-1, 1, 1);
      out(a, b) = out(b, a) = mixed / 24.0;
    }
  }
  return out;
}

} // namespace laplace
} // namespace quadra
