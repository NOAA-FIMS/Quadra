#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

mkdir -p "$BASE/architecture/packages/observation"

cat > "$BASE/architecture/packages/observation/observation_package.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"
#include "../../steps/observation/biomass_index.hpp"
#include "../../steps/observation/catch_age_composition.hpp"
#include "../../../common/model_data.hpp"

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
                  const FleetParameters<T> &parameters,
                  const PopulationState<T> &population,
                  FleetState<T> &fleet) const {
    BiomassIndexPrediction{}(data, parameters, population, fleet);
    CatchAgeCompositionPrediction{}(data, fleet);
  }
};

} // namespace bigeye_v2
CPP

mkdir -p "$BASE/13_observation_package_caa"

cat > "$BASE/13_observation_package_caa/bigeye_v2_13_observation_package_caa_check.cpp" <<'CPP'
#include "../architecture/packages/observation/observation_package.hpp"

#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {
bool nearly_equal(double a, double b, double tol = 1.0e-12) {
  return std::abs(a - b) <= tol;
}
} // namespace

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 2;

  FleetParameters<double> parameters;
  parameters.q_index = 0.01;

  PopulationState<double> population;
  population.spawning_biomass_by_year = {10000.0, 12000.0};

  FleetState<double> fleet;
  fleet.catch_numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    fleet.catch_numbers_at_age[0][a] = static_cast<double>(a + 1);
    fleet.catch_numbers_at_age[1][a] = 2.0 * static_cast<double>(a + 1);
  }

  ObservationPackage{}(data, parameters, population, fleet);

  if (!nearly_equal(fleet.predicted_index_by_year[0], 100.0) ||
      !nearly_equal(fleet.predicted_index_by_year[1], 120.0)) {
    std::cerr << "FAIL: biomass index prediction\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.predicted_catch_age_proportion[0][0], 1.0 / 55.0) ||
      !nearly_equal(fleet.predicted_catch_age_proportion[0][9], 10.0 / 55.0)) {
    std::cerr << std::setprecision(17)
              << "FAIL: age composition prediction got age1 "
              << fleet.predicted_catch_age_proportion[0][0]
              << " age10 "
              << fleet.predicted_catch_age_proportion[0][9]
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA ObservationPackage regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_13_observation_package_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/13_observation_package_caa/bigeye_v2_13_observation_package_caa_check.cpp \
  -o build/examples/bigeye_v2_13_observation_package_caa_check

./build/examples/bigeye_v2_13_observation_package_caa_check
SH

chmod +x run_bigeye_v2_13_observation_package_caa_check.sh

echo "created CAA ObservationPackage"
