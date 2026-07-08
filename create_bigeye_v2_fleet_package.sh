#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

mkdir -p "$BASE/architecture/packages/fleet"

cat > "$BASE/architecture/packages/fleet/fleet_package.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"
#include "../../steps/fleet/baranov_catch.hpp"
#include "../../steps/fleet/fishing_mortality.hpp"
#include "../../steps/fleet/logistic_selectivity.hpp"
#include "../../../common/model_data.hpp"

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
};

} // namespace bigeye_v2
CPP

mkdir -p "$BASE/09_fleet_package_orchestrator_caa"

cat > "$BASE/09_fleet_package_orchestrator_caa/bigeye_v2_09_fleet_package_orchestrator_caa_check.cpp" <<'CPP'
#include "../architecture/packages/fleet/fleet_package.hpp"
#include "../architecture/parameters/fleet_parameters.hpp"
#include "../architecture/parameters/life_history_parameters.hpp"
#include "../architecture/state/fleet_state.hpp"
#include "../architecture/state/life_history_state.hpp"
#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/life_history/life_history.hpp"
#include "../common/model_data.hpp"

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

  LifeHistoryParameters<double> lp;
  lp.log_m_young_offset = std::log(0.75);
  lp.log_m_old_offset = std::log(0.65);

  LifeHistoryState<double> life;
  BigeyeLifeHistory{}(data, lp, life);

  PopulationState<double> pop;
  pop.numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int y = 0; y < data.n_years; ++y) {
    for (int a = 0; a < kAges; ++a) {
      pop.numbers_at_age[static_cast<std::size_t>(y)][a] = 1000.0;
    }
  }

  FleetParameters<double> fp;
  fp.sel_a50 = 5.0;
  fp.sel_slope = 1.0;
  fp.fbar = 0.2;

  FleetState<double> fleet;

  FleetPackage{}(data, fp, life, pop, fleet);

  constexpr double expected_sel_age5 = 0.5;
  constexpr double expected_f_age5 = 0.1;
  constexpr double expected_z_age5 = 0.55;
  constexpr double expected_catch_y0_age1 = 3.0479277436933629;
  constexpr double expected_catch_y0_age5 = 76.918216294456968;
  constexpr double expected_total_catch_y0 = 11918.418702057177;

  if (!nearly_equal(fleet.selectivity_at_age[4], expected_sel_age5)) {
    std::cerr << "FAIL: selectivity age5\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.f_at_age[4], expected_f_age5)) {
    std::cerr << "FAIL: F age5\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.z_at_age[4], expected_z_age5)) {
    std::cerr << "FAIL: Z age5\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.catch_numbers_at_age[0][0],
                    expected_catch_y0_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch y0 age1 got "
              << fleet.catch_numbers_at_age[0][0]
              << " expected " << expected_catch_y0_age1
              << " diff "
              << (fleet.catch_numbers_at_age[0][0] -
                  expected_catch_y0_age1)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.catch_numbers_at_age[0][4],
                    expected_catch_y0_age5)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch y0 age5 got "
              << fleet.catch_numbers_at_age[0][4]
              << " expected " << expected_catch_y0_age5
              << " diff "
              << (fleet.catch_numbers_at_age[0][4] -
                  expected_catch_y0_age5)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.total_catch_biomass_by_year[0],
                    expected_total_catch_y0)) {
    std::cerr << std::setprecision(17)
              << "FAIL: total catch y0 got "
              << fleet.total_catch_biomass_by_year[0]
              << " expected " << expected_total_catch_y0
              << " diff "
              << (fleet.total_catch_biomass_by_year[0] -
                  expected_total_catch_y0)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA FleetPackage regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_09_fleet_package_orchestrator_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/09_fleet_package_orchestrator_caa/bigeye_v2_09_fleet_package_orchestrator_caa_check.cpp \
  -o build/examples/bigeye_v2_09_fleet_package_orchestrator_caa_check

./build/examples/bigeye_v2_09_fleet_package_orchestrator_caa_check
SH

chmod +x run_bigeye_v2_09_fleet_package_orchestrator_caa_check.sh

echo "created CAA FleetPackage"
