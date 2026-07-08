#pragma once

#include "../../parameters/population_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"

namespace bigeye_v2 {

template <typename T>
struct PopulationContext {
  const PopulationParameters<T> *parameters = nullptr;
  const LifeHistoryState<T> *life = nullptr;
  const FleetState<T> *fleet = nullptr;
  PopulationState<T> *population = nullptr;
};

} // namespace bigeye_v2
