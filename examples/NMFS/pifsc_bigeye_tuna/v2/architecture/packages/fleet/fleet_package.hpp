#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"
#include "../../steps/fleet/baranov_catch.hpp"
#include "../../steps/fleet/fishing_mortality.hpp"
#include "../../steps/fleet/logistic_selectivity.hpp"
#include "fleet_context.hpp"

namespace bigeye_v2 {

//------------------------------------------------------------
// FleetPackage
//
// Purpose
// -------
// Composes the fleet steps for one fleet.
//
// Sequence
// --------
// Select
// Fish
// Catch
//
// Consumes
// --------
// BigeyeModelData
// FleetParameters
// LifeHistoryState
// PopulationState
//
// Produces
// --------
// FleetState::selectivity_at_age
// FleetState::f_at_age
// FleetState::z_at_age
// FleetState::catch_numbers_at_age
// FleetState::catch_biomass_at_age
// FleetState::total_catch_biomass_by_year
//
// Notes
// -----
// State owns memory.
// Steps own algorithms.
// Packages orchestrate related steps.
//------------------------------------------------------------
struct FleetPackage {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const LifeHistoryState<T> &life,
                  const PopulationState<T> &population,
                  FleetState<T> &fleet) const {
    LogisticSelectivity{}(data, parameters, fleet);
    FishingMortality{}(parameters, life, fleet);
    BaranovCatch{}(data, life, population, fleet);
  }

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetContext<T> &context) const {
    LogisticSelectivity{}(data, *context.parameters, *context.fleet);
    FishingMortality{}(*context.parameters, *context.life, *context.fleet);
    BaranovCatch{}(data, *context.life, *context.population, *context.fleet);
  }
};

} // namespace bigeye_v2
