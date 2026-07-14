#pragma once

#include "../../state/population_state.hpp"

namespace bigeye_v2 {

// Accumulates survivors into the terminal plus group.
struct PlusGroup {
  template <typename T> void operator()(PopulationState<T> &population) const {
    const auto ny = population.survivors_at_age.size();

    for (std::size_t y = 0; y + 1 < ny; ++y) {
      constexpr int plus = kAges - 1;

      population.numbers_at_age[y + 1][plus] =
          population.survivors_at_age[y][plus - 1] +
          population.survivors_at_age[y][plus];
    }
  }
};

} // namespace bigeye_v2
