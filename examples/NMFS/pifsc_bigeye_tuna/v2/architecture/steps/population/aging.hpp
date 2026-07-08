#pragma once

#include "../../state/population_state.hpp"

#include <array>

namespace bigeye_v2 {

//------------------------------------------------------------
// Aging
//
// Purpose
// -------
// Advances survivors into the next age class for the next year.
//
// Consumes
// --------
// PopulationState::survivors_at_age
//
// Produces
// --------
// PopulationState::numbers_at_age for year + 1
//
// Notes
// -----
// Stateless.
// Owns no memory beyond resizing population-owned state.
// Does not apply recruitment.
// Does not handle the plus group.
//------------------------------------------------------------
struct Aging {
  template <typename T>
  void operator()(PopulationState<T> &population) const {
    const auto ny = population.survivors_at_age.size();

    if (population.numbers_at_age.size() < ny) {
      population.numbers_at_age.assign(ny, std::array<T, kAges>{});
    }

    for (std::size_t y = 0; y + 1 < ny; ++y) {
      for (int a = 1; a < kAges - 1; ++a) {
        population.numbers_at_age[y + 1][a] =
            population.survivors_at_age[y][a - 1];
      }
    }
  }
};

} // namespace bigeye_v2
