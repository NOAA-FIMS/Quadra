#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"

#include <vector>

namespace bigeye_v2 {

// Predicts a biomass index from population spawning biomass and fleet q.
struct BiomassIndexPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const PopulationState<T> &population,
                  FleetState<T> &fleet) const {
    fleet.predicted_index_by_year.assign(static_cast<std::size_t>(data.n_years),
                                         T(0.0));

    for (std::size_t y = 0; y < population.spawning_biomass_by_year.size();
         ++y) {
      fleet.predicted_index_by_year[y] =
          parameters.q_index * population.spawning_biomass_by_year[y];
    }
  }
};

} // namespace bigeye_v2
