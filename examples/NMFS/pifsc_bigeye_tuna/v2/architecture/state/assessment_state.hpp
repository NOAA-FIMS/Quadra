#pragma once

#include "fleet_state.hpp"
#include "life_history_state.hpp"
#include "likelihood_state.hpp"
#include "population_state.hpp"

#include <vector>

namespace bigeye_v2 {

template <typename T>
struct AssessmentState {
  LifeHistoryState<T> life;
  std::vector<PopulationState<T>> populations;
  std::vector<FleetState<T>> fleets;
  LikelihoodState<T> likelihood;
};

} // namespace bigeye_v2
