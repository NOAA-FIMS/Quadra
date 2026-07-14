#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"
#include "../../steps/observation/biomass_index.hpp"
#include "../../steps/observation/catch_age_composition.hpp"
#include "observation_context.hpp"

namespace bigeye_v2 {

//------------------------------------------------------------
// ObservationPackage
//
// Purpose
// -------
// Composes observation prediction steps for one fleet.
//
// Sequence
// --------
// Predict biomass index
// Predict catch age composition
//
// Notes
// -----
// State owns memory.
// Steps own algorithms.
// Packages orchestrate related steps.
//------------------------------------------------------------
struct ObservationPackage {

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const ObservationContext<T> &context) const {
    BiomassIndexPrediction{}(data, *context.parameters, *context.population,
                             *context.fleet);
    CatchAgeCompositionPrediction{}(data, *context.fleet);
  }

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const PopulationState<T> &population,
                  FleetState<T> &fleet) const {
    BiomassIndexPrediction{}(data, parameters, population, fleet);
    CatchAgeCompositionPrediction{}(data, fleet);
  }
};

} // namespace bigeye_v2
