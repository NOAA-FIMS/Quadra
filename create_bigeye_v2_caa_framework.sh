#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2/architecture"

mkdir -p "$BASE/data" \
         "$BASE/parameters" \
         "$BASE/state" \
         "$BASE/steps/life_history" \
         "$BASE/steps/population" \
         "$BASE/steps/movement" \
         "$BASE/steps/fleet" \
         "$BASE/steps/observation" \
         "$BASE/steps/likelihood" \
         "$BASE/assessment"

cat > "$BASE/state/life_history_state.hpp" <<'CPP'
#pragma once

#include "../../common/bigeye_constants.hpp"

#include <array>

namespace bigeye_v2 {

template <typename T>
struct LifeHistoryState {
  std::array<T, kAges> m_at_age{};
  std::array<T, kAges> weight_at_age{};
  std::array<T, kAges> maturity_at_age{};
};

} // namespace bigeye_v2
CPP

cat > "$BASE/state/population_state.hpp" <<'CPP'
#pragma once

#include "../../common/bigeye_constants.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T>
struct PopulationState {
  std::vector<T> recruits_by_year{};
  std::vector<std::array<T, kAges>> numbers_at_age{};
  std::vector<T> spawning_biomass_by_year{};
};

} // namespace bigeye_v2
CPP

cat > "$BASE/state/fleet_state.hpp" <<'CPP'
#pragma once

#include "../../common/bigeye_constants.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T>
struct FleetState {
  std::array<T, kAges> selectivity_at_age{};
  std::array<T, kAges> f_at_age{};
  std::array<T, kAges> z_at_age{};

  std::vector<std::array<T, kAges>> catch_numbers_at_age{};
  std::vector<std::array<T, kAges>> catch_biomass_at_age{};
  std::vector<T> total_catch_biomass_by_year{};
};

} // namespace bigeye_v2
CPP

cat > "$BASE/state/likelihood_state.hpp" <<'CPP'
#pragma once

namespace bigeye_v2 {

template <typename T>
struct LikelihoodState {
  T catch_nll = T(0.0);
  T index_nll = T(0.0);
  T agecomp_nll = T(0.0);
  T total_nll = T(0.0);
};

} // namespace bigeye_v2
CPP

cat > "$BASE/state/assessment_state.hpp" <<'CPP'
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
CPP

cat > "$BASE/parameters/life_history_parameters.hpp" <<'CPP'
#pragma once

namespace bigeye_v2 {

template <typename T>
struct LifeHistoryParameters {
  T log_m_young_offset = T(0.0);
  T log_m_old_offset = T(0.0);
};

} // namespace bigeye_v2
CPP

cat > "$BASE/parameters/population_parameters.hpp" <<'CPP'
#pragma once

namespace bigeye_v2 {

template <typename T>
struct PopulationParameters {
  T r0 = T(1000.0);
};

} // namespace bigeye_v2
CPP

cat > "$BASE/parameters/fleet_parameters.hpp" <<'CPP'
#pragma once

namespace bigeye_v2 {

template <typename T>
struct FleetParameters {
  T sel_a50 = T(5.0);
  T sel_slope = T(1.0);
  T fbar = T(0.2);
  T q_index = T(0.01);
  T catch_sigma = T(0.1);
  T index_sigma = T(0.2);
};

} // namespace bigeye_v2
CPP

cat > "$BASE/parameters/movement_parameters.hpp" <<'CPP'
#pragma once

namespace bigeye_v2 {

template <typename T>
struct MovementParameters {
  T move_0_to_1 = T(0.0);
  T move_1_to_0 = T(0.0);
};

} // namespace bigeye_v2
CPP

cat > "$BASE/parameters/assessment_parameters.hpp" <<'CPP'
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
CPP

cat > "$BASE/assessment/assessment_cycle.hpp" <<'CPP'
#pragma once

namespace bigeye_v2 {

// The assessment cycle owns execution order.
// Scientific steps own equations.
struct AssessmentCycle {
};

} // namespace bigeye_v2
CPP

echo "created CAA framework scaffold"
