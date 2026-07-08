#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level12_timestep_check"

cat > "$BASE/common/timestep.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"
#include "fleet.hpp"
#include "model_data.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace bigeye_v2 {

template <typename T>
struct TimeStepPopulationState {
  std::vector<std::array<T, kAges>> numbers_at_age{};
};

struct AnnualTimeStep {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeDerived<T> &life,
                  const BigeyeFleetDerived<T> &fleet,
                  const T &recruitment,
                  TimeStepPopulationState<T> &pop) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    if (pop.numbers_at_age.size() < ny) {
      pop.numbers_at_age.assign(ny, std::array<T, kAges>{});
    }

    for (std::size_t y = 0; y + 1 < ny; ++y) {
      pop.numbers_at_age[y + 1][0] = recruitment;

      for (int a = 1; a < kAges - 1; ++a) {
        pop.numbers_at_age[y + 1][a] =
            pop.numbers_at_age[y][a - 1] *
            std::exp(-fleet.z_at_age[a - 1]);
      }

      const int plus = kAges - 1;

      pop.numbers_at_age[y + 1][plus] =
          pop.numbers_at_age[y][plus - 1] *
              std::exp(-fleet.z_at_age[plus - 1]) +
          pop.numbers_at_age[y][plus] *
              std::exp(-fleet.z_at_age[plus]);
    }
  }
};

struct AnnualCatchAtStartOfYear {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeDerived<T> &life,
                  const BigeyeFleetDerived<T> &fleet_state,
                  const TimeStepPopulationState<T> &pop,
                  BigeyeFleetDerived<T> &fleet) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    fleet.catch_numbers_at_age.assign(ny, std::array<T, kAges>{});
    fleet.catch_biomass_at_age.assign(ny, std::array<T, kAges>{});
    fleet.total_catch_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        const T z = fleet_state.z_at_age[a];

        const T cn =
            pop.numbers_at_age[y][a] *
            (fleet_state.f_at_age[a] / z) *
            (T(1.0) - std::exp(-z));

        fleet.catch_numbers_at_age[y][a] = cn;
        fleet.catch_biomass_at_age[y][a] = cn * life.weight_at_age[a];
        fleet.total_catch_biomass_by_year[y] +=
            fleet.catch_biomass_at_age[y][a];
      }
    }
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/level12_timestep_check/bigeye_v2_level12_timestep_check.cpp" <<'CPP'
#include "../common/fleet.hpp"
#include "../common/life_history.hpp"
#include "../common/timestep.hpp"

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

  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);

  BigeyeDerived<double> life;
  BigeyeLifeHistory{}(data, p, life);

  BigeyeFleetParameters<double> fp;
  fp.sel_a50 = 5.0;
  fp.sel_slope = 1.0;
  fp.fbar = 0.2;

  BigeyeFleetDerived<double> fleet;
  FleetLogisticSelectivity{}(data, fp, fleet);
  FleetFishingMortality{}(data, fp, life, fleet);

  TimeStepPopulationState<double> pop;
  pop.numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  // Initial year numbers-at-age.
  for (int a = 0; a < kAges; ++a) {
    pop.numbers_at_age[0][a] = 1000.0;
  }

  AnnualTimeStep{}(data, life, fleet, 1000.0, pop);
  AnnualCatchAtStartOfYear{}(data, life, fleet, pop, fleet);

  constexpr double expected_y1_age1 = 1000.0;
  constexpr double expected_y1_age2 = 711.9177620339338;
  constexpr double expected_y1_age10 = 1226.9201954177395;
  constexpr double expected_catch_y0_age1 = 3.0479277436933629;
  constexpr double expected_catch_y1_age2 = 6.89337826858519;

  if (!nearly_equal(pop.numbers_at_age[1][0], expected_y1_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age1 got " << pop.numbers_at_age[1][0]
              << " expected " << expected_y1_age1 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.numbers_at_age[1][1], expected_y1_age2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age2 got " << pop.numbers_at_age[1][1]
              << " expected " << expected_y1_age2
              << " diff " << (pop.numbers_at_age[1][1] - expected_y1_age2)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.numbers_at_age[1][9], expected_y1_age10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age10 got " << pop.numbers_at_age[1][9]
              << " expected " << expected_y1_age10
              << " diff " << (pop.numbers_at_age[1][9] - expected_y1_age10)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.catch_numbers_at_age[0][0],
                    expected_catch_y0_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch y0 age1 got "
              << fleet.catch_numbers_at_age[0][0]
              << " expected " << expected_catch_y0_age1
              << " diff "
              << (fleet.catch_numbers_at_age[0][0] - expected_catch_y0_age1)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.catch_numbers_at_age[1][1],
                    expected_catch_y1_age2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch y1 age2 got "
              << fleet.catch_numbers_at_age[1][1]
              << " expected " << expected_catch_y1_age2
              << " diff "
              << (fleet.catch_numbers_at_age[1][1] - expected_catch_y1_age2)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level12 annual timestep regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level12_timestep_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level12_timestep_check/bigeye_v2_level12_timestep_check.cpp \
  -o build/examples/bigeye_v2_level12_timestep_check

./build/examples/bigeye_v2_level12_timestep_check
SH

chmod +x run_bigeye_v2_level12_timestep_check.sh

echo "created Bigeye v2 Level12 annual timestep check"
