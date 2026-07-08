#pragma once

#include "bigeye_constants.hpp"
#include "model_data.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T>
struct BigeyePopulationDerived {
  std::vector<std::array<T, kAges>> numbers_at_age{};
};

template <typename T>
struct BigeyeMovementParameters {
  T move_0_to_1 = T(0.10);
  T move_1_to_0 = T(0.05);
};

struct TwoPopulationMovement {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeMovementParameters<T> &p,
                  std::vector<BigeyePopulationDerived<T>> &pop) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    if (pop.size() != 2) {
      return;
    }

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        const T n0 = pop[0].numbers_at_age[y][a];
        const T n1 = pop[1].numbers_at_age[y][a];

        pop[0].numbers_at_age[y][a] =
            n0 * (T(1.0) - p.move_0_to_1) + n1 * p.move_1_to_0;

        pop[1].numbers_at_age[y][a] =
            n1 * (T(1.0) - p.move_1_to_0) + n0 * p.move_0_to_1;
      }
    }
  }
};

} // namespace bigeye_v2
