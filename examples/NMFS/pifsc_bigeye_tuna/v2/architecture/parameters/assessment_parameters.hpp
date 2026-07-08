#pragma once

#include "fleet_parameters.hpp"
#include "life_history_parameters.hpp"
#include "movement_parameters.hpp"
#include "population_parameters.hpp"

#include <vector>

namespace bigeye_v2 {

template <typename T>
struct AssessmentParameters {
  LifeHistoryParameters<T> life;
  std::vector<PopulationParameters<T>> populations;
  std::vector<FleetParameters<T>> fleets;
  MovementParameters<T> movement;
};

} // namespace bigeye_v2
