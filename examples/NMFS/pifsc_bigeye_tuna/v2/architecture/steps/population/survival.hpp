#pragma once

#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace bigeye_v2 {

//------------------------------------------------------------
// Survival
//
// Purpose
// -------
// Computes survivors-at-age after total mortality.
//
// Consumes
// --------
// PopulationState::numbers_at_age
// FleetState::z_at_age
//
// Produces
// --------
// PopulationState::survivors_at_age
//
// Notes
// -----
// Stateless.
// Owns no memory beyond resizing population-owned state.
// Does not age fish.
// Does not apply recruitment.
// Does not handle the plus group.
//------------------------------------------------------------
struct Survival {
  template <typename T>
  void operator()(const FleetState<T> &fleet,
                  PopulationState<T> &population) const {
    const auto ny = population.numbers_at_age.size();

    population.survivors_at_age.assign(ny, std::array<T, kAges>{});

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        population.survivors_at_age[y][a] =
            population.numbers_at_age[y][a] * std::exp(-fleet.z_at_age[a]);
      }
    }
  }
};

} // namespace bigeye_v2
