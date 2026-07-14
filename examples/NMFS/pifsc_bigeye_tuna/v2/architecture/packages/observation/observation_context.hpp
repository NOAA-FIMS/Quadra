#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"

namespace bigeye_v2 {

template <typename T> struct ObservationContext {
  const FleetParameters<T> *parameters = nullptr;
  const PopulationState<T> *population = nullptr;
  FleetState<T> *fleet = nullptr;
};

} // namespace bigeye_v2
