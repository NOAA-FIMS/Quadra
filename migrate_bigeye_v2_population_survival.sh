#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

cat > "$BASE/architecture/steps/population/survival.hpp" <<'CPP'
#pragma once

#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace bigeye_v2 {

//------------------------------------------------------------
// Survival
//
// Purpose
// -------
// Computes survivors-at-age after total mortality.
//
// Consumes
// --------
// PopulationState::numbers_at_age
// FleetState::z_at_age
//
// Produces
// --------
// PopulationState::survivors_at_age
//
// Notes
// -----
// Stateless.
// Owns no memory beyond resizing population-owned state.
// Does not age fish.
// Does not apply recruitment.
// Does not handle the plus group.
//------------------------------------------------------------
struct Survival {
  template <typename T>
  void operator()(const FleetState<T> &fleet,
                  PopulationState<T> &population) const {
    const auto ny = population.numbers_at_age.size();

    population.survivors_at_age.assign(ny, std::array<T, kAges>{});

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        population.survivors_at_age[y][a] =
            population.numbers_at_age[y][a] * std::exp(-fleet.z_at_age[a]);
      }
    }
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/state/population_state.hpp")
s = p.read_text()
if "survivors_at_age" not in s:
    s = s.replace(
        "  std::vector<std::array<T, kAges>> numbers_at_age{};",
        "  std::vector<std::array<T, kAges>> numbers_at_age{};\n"
        "  std::vector<std::array<T, kAges>> survivors_at_age{};"
    )
p.write_text(s)
PY

mkdir -p "$BASE/05_population_survival_caa"

cat > "$BASE/05_population_survival_caa/bigeye_v2_05_population_survival_caa_check.cpp" <<'CPP'
#include "../architecture/state/fleet_state.hpp"
#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/population/survival.hpp"

#include <array>
#include <cmath>
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

  PopulationState<double> population;
  population.numbers_at_age.assign(2, std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    population.numbers_at_age[0][a] = 1000.0;
    population.numbers_at_age[1][a] = 2000.0;
  }

  FleetState<double> fleet;
  for (int a = 0; a < kAges; ++a) {
    fleet.z_at_age[a] = 0.1 * static_cast<double>(a + 1);
  }

  Survival{}(fleet, population);

  constexpr double expected_y0_age1 = 904.83741803595957;
  constexpr double expected_y0_age10 = 367.87944117144235;
  constexpr double expected_y1_age1 = 1809.6748360719191;

  if (!nearly_equal(population.survivors_at_age[0][0], expected_y0_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y0 age1 survivor got "
              << population.survivors_at_age[0][0]
              << " expected " << expected_y0_age1
              << " diff " << (population.survivors_at_age[0][0] -
                               expected_y0_age1)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.survivors_at_age[0][9], expected_y0_age10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y0 age10 survivor got "
              << population.survivors_at_age[0][9]
              << " expected " << expected_y0_age10
              << " diff " << (population.survivors_at_age[0][9] -
                               expected_y0_age10)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.survivors_at_age[1][0], expected_y1_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age1 survivor got "
              << population.survivors_at_age[1][0]
              << " expected " << expected_y1_age1
              << " diff " << (population.survivors_at_age[1][0] -
                               expected_y1_age1)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA population survival regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_05_population_survival_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/05_population_survival_caa/bigeye_v2_05_population_survival_caa_check.cpp \
  -o build/examples/bigeye_v2_05_population_survival_caa_check

./build/examples/bigeye_v2_05_population_survival_caa_check
SH

chmod +x run_bigeye_v2_05_population_survival_caa_check.sh

echo "migrated CAA population survival step"
