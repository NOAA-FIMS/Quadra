#pragma once

#include "../packages/fleet/fleet_context.hpp"
#include "../packages/fleet/fleet_package.hpp"
#include "../packages/life_history/life_history_package.hpp"
#include "../packages/likelihood/likelihood_context.hpp"
#include "../packages/likelihood/likelihood_package.hpp"
#include "../packages/observation/observation_context.hpp"
#include "../packages/observation/observation_package.hpp"
#include "../packages/population/population_context.hpp"
#include "../packages/population/population_package.hpp"
#include "../parameters/assessment_parameters.hpp"
#include "../state/assessment_state.hpp"
#include "../../common/model_data.hpp"

namespace bigeye_v2 {

// AssessmentCycle owns orchestration.
// Packages own workflows.
// Steps own algorithms.
// State owns memory.
struct AssessmentCycle {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const AssessmentParameters<T> &parameters,
                  AssessmentState<T> &state) const {
    LifeHistoryPackage{}(data, parameters.life, state.life);

    if (state.populations.size() < parameters.populations.size()) {
      state.populations.resize(parameters.populations.size());
    }

    if (state.fleets.size() < parameters.fleets.size()) {
      state.fleets.resize(parameters.fleets.size());
    }

    if (parameters.populations.empty() || parameters.fleets.empty()) {
      return;
    }

    FleetContext<T> fleet_context{
        &parameters.fleets[0],
        &state.life,
        &state.populations[0],
        &state.fleets[0]};

    PopulationContext<T> population_context{
        &parameters.populations[0],
        &state.life,
        &state.fleets[0],
        &state.populations[0]};

    ObservationContext<T> observation_context{
        &parameters.fleets[0],
        &state.populations[0],
        &state.fleets[0]};

    LikelihoodContext<T> likelihood_context{
        &parameters.fleets[0],
        &state.fleets[0],
        &state.likelihood};

    // First pass initializes fleet mortality from current population state.
    FleetPackage{}(data, fleet_context);

    // Population advances using fleet total mortality.
    PopulationPackage{}(data, population_context);

    // Recompute fleet predictions from advanced population state.
    FleetPackage{}(data, fleet_context);

    // Prediction layer.
    ObservationPackage{}(data, observation_context);

    // Objective layer.
    LikelihoodPackage{}(data, likelihood_context);
  }
};

} // namespace bigeye_v2
