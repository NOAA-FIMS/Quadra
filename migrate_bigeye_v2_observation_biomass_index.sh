#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

mkdir -p "$BASE/architecture/steps/observation"

cat > "$BASE/architecture/steps/observation/biomass_index.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"
#include "../../../common/model_data.hpp"

#include <vector>

namespace bigeye_v2 {

// Predicts a biomass index from population spawning biomass and fleet q.
struct BiomassIndexPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const PopulationState<T> &population,
                  FleetState<T> &fleet) const {
    fleet.predicted_index_by_year.assign(
        static_cast<std::size_t>(data.n_years), T(0.0));

    for (std::size_t y = 0; y < population.spawning_biomass_by_year.size(); ++y) {
      fleet.predicted_index_by_year[y] =
          parameters.q_index * population.spawning_biomass_by_year[y];
    }
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path
p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/state/fleet_state.hpp")
s = p.read_text()
if "predicted_index_by_year" not in s:
    s = s.replace(
        "  std::vector<T> total_catch_biomass_by_year{};",
        "  std::vector<T> total_catch_biomass_by_year{};\n"
        "  std::vector<T> predicted_index_by_year{};"
    )
p.write_text(s)
PY

mkdir -p "$BASE/11_observation_biomass_index_caa"

cat > "$BASE/11_observation_biomass_index_caa/bigeye_v2_11_observation_biomass_index_caa_check.cpp" <<'CPP'
#include "../architecture/parameters/fleet_parameters.hpp"
#include "../architecture/state/fleet_state.hpp"
#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/observation/biomass_index.hpp"
#include "../common/model_data.hpp"

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
  data.n_years = 3;

  PopulationState<double> population;
  population.spawning_biomass_by_year = {10000.0, 12000.0, 9000.0};

  FleetParameters<double> parameters;
  parameters.q_index = 0.01;

  FleetState<double> fleet;

  BiomassIndexPrediction{}(data, parameters, population, fleet);

  constexpr double expected_y0 = 100.0;
  constexpr double expected_y1 = 120.0;
  constexpr double expected_y2 = 90.0;

  if (!nearly_equal(fleet.predicted_index_by_year[0], expected_y0) ||
      !nearly_equal(fleet.predicted_index_by_year[1], expected_y1) ||
      !nearly_equal(fleet.predicted_index_by_year[2], expected_y2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: predicted index values got "
              << fleet.predicted_index_by_year[0] << ", "
              << fleet.predicted_index_by_year[1] << ", "
              << fleet.predicted_index_by_year[2] << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA biomass index observation regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_11_observation_biomass_index_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/11_observation_biomass_index_caa/bigeye_v2_11_observation_biomass_index_caa_check.cpp \
  -o build/examples/bigeye_v2_11_observation_biomass_index_caa_check

./build/examples/bigeye_v2_11_observation_biomass_index_caa_check
SH

chmod +x run_bigeye_v2_11_observation_biomass_index_caa_check.sh

echo "migrated CAA biomass index observation step"
