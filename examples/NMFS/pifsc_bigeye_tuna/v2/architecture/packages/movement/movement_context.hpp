#pragma once

#include "../../parameters/movement_parameters.hpp"
#include "../../state/population_state.hpp"

#include <vector>

namespace bigeye_v2 {

template <typename T>
struct MovementContext {
  const MovementParameters<T> *parameters = nullptr;
  std::vector<PopulationState<T>> *populations = nullptr;
};

} // namespace bigeye_v2
