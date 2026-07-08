#!/usr/bin/env bash
set -euo pipefail

p="examples/NMFS/pifsc_bigeye_tuna/v2/architecture/assessment/assessment_cycle.hpp"

cat > "$p" <<'CPP'
#pragma once

#include "../packages/fleet/fleet_context.hpp"
#include "../packages/fleet/fleet_package.hpp"
#include "../packages/life_history/life_history_package.hpp"
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

    // First pass initializes fleet mortality from current population state.
    FleetPackage{}(data, fleet_context);

    // Population advances using fleet total mortality.
    PopulationPackage{}(data, population_context);

    // Recompute fleet predictions from advanced population state.
    FleetPackage{}(data, fleet_context);

    // Prediction-only observation layer.
    ObservationPackage{}(data, observation_context);
  }
};

} // namespace bigeye_v2
CPP

echo "updated AssessmentCycle to use package contexts"
