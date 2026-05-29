#include <Eigen/Sparse>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/hdot_active_direction_discovery.hpp"

int main() {
  std::vector<Eigen::SparseMatrix<double>> Hdots(4);

  Hdots[0].resize(3, 3);
  Hdots[0].makeCompressed();

  {
    std::vector<Eigen::Triplet<double>> t;
    t.emplace_back(0, 0, 1.0);
    t.emplace_back(1, 1, 2.0);
    Hdots[1].resize(3, 3);
    Hdots[1].setFromTriplets(t.begin(), t.end());
  }

  {
    std::vector<Eigen::Triplet<double>> t;
    t.emplace_back(0, 1, 1.0e-8);
    Hdots[2].resize(3, 3);
    Hdots[2].setFromTriplets(t.begin(), t.end());
  }

  {
    std::vector<Eigen::Triplet<double>> t;
    t.emplace_back(2, 2, 5.0);
    Hdots[3].resize(3, 3);
    Hdots[3].setFromTriplets(t.begin(), t.end());
  }

  const auto result =
      quadra::laplace::discover_active_directions_from_hdots(Hdots, 1.0e-6);

  if (result.active_directions.size() != 2 ||
      result.active_directions[0] != 1 || result.active_directions[1] != 3) {
    throw std::runtime_error("active direction discovery returned wrong set.");
  }

  if (result.hdot_norms.size() != 4) {
    throw std::runtime_error("wrong norm vector size.");
  }

  std::cout << "Hdot active-direction discovery tests passed\n";
  return 0;
}
