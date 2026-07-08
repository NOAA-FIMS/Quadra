#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

mkdir -p "$BASE/architecture/packages/life_history"
mkdir -p "$BASE/architecture/assessment"

cat > "$BASE/architecture/packages/life_history/life_history_package.hpp" <<'CPP'
#pragma once

#include "../../parameters/life_history_parameters.hpp"
#include "../../state/life_history_state.hpp"
#include "../../steps/life_history/life_history.hpp"
#include "../../../common/model_data.hpp"

namespace bigeye_v2 {

struct LifeHistoryPackage {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const LifeHistoryParameters<T> &parameters,
                  LifeHistoryState<T> &life) const {
    BigeyeLifeHistory{}(data, parameters, life);
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/architecture/assessment/assessment_cycle.hpp" <<'CPP'
#pragma once

#include "../packages/fleet/fleet_package.hpp"
#include "../packages/life_history/life_history_package.hpp"
#include "../packages/population/population_package.hpp"
#include "../parameters/assessment_parameters.hpp"
#include "../state/assessment_state.hpp"
#include "../../common/model_data.hpp"

namespace bigeye_v2 {

// AssessmentCycle orchestrates packages.
// Packages orchestrate steps.
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

    // Minimal first pass:
    // one population, one fleet, no movement yet.
    //
    // The population package needs total mortality. Until the movement /
    // multi-fleet mortality aggregator exists, compute fleet state first
    // from any initialized population state, then use it to advance population.
    if (!parameters.populations.empty() && !parameters.fleets.empty()) {
      FleetPackage{}(
          data,
          parameters.fleets[0],
          state.life,
          state.populations[0],
          state.fleets[0]);

      PopulationPackage{}(
          data,
          parameters.populations[0],
          state.life,
          state.fleets[0],
          state.populations[0]);

      FleetPackage{}(
          data,
          parameters.fleets[0],
          state.life,
          state.populations[0],
          state.fleets[0]);
    }
  }
};

} // namespace bigeye_v2
CPP

mkdir -p "$BASE/10_assessment_cycle_caa"

cat > "$BASE/10_assessment_cycle_caa/bigeye_v2_10_assessment_cycle_caa_check.cpp" <<'CPP'
#include "../architecture/assessment/assessment_cycle.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {
bool nearly_equal(double a, double b, double tol = 1.0e-6) {
  return std::abs(a - b) <= tol;
}
} // namespace

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 3;

  AssessmentParameters<double> parameters;
  parameters.life.log_m_young_offset = std::log(0.75);
  parameters.life.log_m_old_offset = std::log(0.65);

  parameters.populations.resize(1);
  parameters.populations[0].r0 = 1000.0;

  parameters.fleets.resize(1);
  parameters.fleets[0].sel_a50 = 5.0;
  parameters.fleets[0].sel_slope = 1.0;
  parameters.fleets[0].fbar = 0.2;

  AssessmentState<double> state;
  state.populations.resize(1);
  state.fleets.resize(1);

  state.populations[0].numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    state.populations[0].numbers_at_age[0][a] = 1000.0;
  }

  AssessmentCycle{}(data, parameters, state);

  constexpr double expected_m_age1 = 0.3375;
  constexpr double expected_pop_y1_age2 = 710.98976678322833;
  constexpr double expected_fleet_z_age5 = 0.55;
  constexpr double expected_total_catch_y0 = 11918.418702057177;

  if (!nearly_equal(state.life.m_at_age[0], expected_m_age1)) {
    std::cerr << "FAIL: life history not computed\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.populations[0].numbers_at_age[1][1],
                    expected_pop_y1_age2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: population y1 age2 got "
              << state.populations[0].numbers_at_age[1][1]
              << " expected " << expected_pop_y1_age2
              << " diff "
              << (state.populations[0].numbers_at_age[1][1] -
                  expected_pop_y1_age2)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.fleets[0].z_at_age[4], expected_fleet_z_age5)) {
    std::cerr << "FAIL: fleet mortality not computed\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.fleets[0].total_catch_biomass_by_year[0],
                    expected_total_catch_y0)) {
    std::cerr << std::setprecision(17)
              << "FAIL: total catch y0 got "
              << state.fleets[0].total_catch_biomass_by_year[0]
              << " expected " << expected_total_catch_y0
              << " diff "
              << (state.fleets[0].total_catch_biomass_by_year[0] -
                  expected_total_catch_y0)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA AssessmentCycle regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_10_assessment_cycle_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/10_assessment_cycle_caa/bigeye_v2_10_assessment_cycle_caa_check.cpp \
  -o build/examples/bigeye_v2_10_assessment_cycle_caa_check

./build/examples/bigeye_v2_10_assessment_cycle_caa_check
SH

chmod +x run_bigeye_v2_10_assessment_cycle_caa_check.sh

echo "created CAA AssessmentCycle"
