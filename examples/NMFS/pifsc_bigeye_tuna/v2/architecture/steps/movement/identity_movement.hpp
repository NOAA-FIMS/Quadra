#pragma once

#include "../../parameters/movement_parameters.hpp"
#include "../../state/population_state.hpp"

#include <vector>

namespace bigeye_v2 {

// Identity movement.
// Leaves population numbers unchanged.
struct IdentityMovement {
  template <typename T>
  void operator()(const MovementParameters<T> &,
                  std::vector<PopulationState<T>> &) const {
    // Intentionally empty.
  }
};

} // namespace bigeye_v2
