#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

cat > "$BASE/architecture/steps/observation/catch_age_composition.hpp" <<'CPP'
#pragma once

#include "../../state/fleet_state.hpp"
#include "../../../common/bigeye_constants.hpp"
#include "../../../common/model_data.hpp"

#include <array>

namespace bigeye_v2 {

// Predicts catch age composition from fleet catch numbers-at-age.
struct CatchAgeCompositionPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  FleetState<T> &fleet) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    fleet.predicted_catch_age_proportion.assign(
        ny, std::array<T, kAges>{});

    for (std::size_t y = 0; y < ny; ++y) {
      T total = T(0.0);

      for (int a = 0; a < kAges; ++a) {
        total += fleet.catch_numbers_at_age[y][a];
      }

      for (int a = 0; a < kAges; ++a) {
        fleet.predicted_catch_age_proportion[y][a] =
            fleet.catch_numbers_at_age[y][a] / total;
      }
    }
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path
p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/state/fleet_state.hpp")
s = p.read_text()
if "predicted_catch_age_proportion" not in s:
    s = s.replace(
        "  std::vector<T> predicted_index_by_year{};",
        "  std::vector<T> predicted_index_by_year{};\n"
        "  std::vector<std::array<T, kAges>> predicted_catch_age_proportion{};"
    )
p.write_text(s)
PY

mkdir -p "$BASE/12_observation_agecomp_caa"

cat > "$BASE/12_observation_agecomp_caa/bigeye_v2_12_observation_agecomp_caa_check.cpp" <<'CPP'
#include "../architecture/state/fleet_state.hpp"
#include "../architecture/steps/observation/catch_age_composition.hpp"
#include "../common/model_data.hpp"

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

  FleetState<double> fleet;
  fleet.catch_numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    fleet.catch_numbers_at_age[0][a] = static_cast<double>(a + 1);
    fleet.catch_numbers_at_age[1][a] = 2.0 * static_cast<double>(a + 1);
  }

  CatchAgeCompositionPrediction{}(data, fleet);

  constexpr double expected_age1 = 1.0 / 55.0;
  constexpr double expected_age10 = 10.0 / 55.0;

  if (!nearly_equal(fleet.predicted_catch_age_proportion[0][0],
                    expected_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: age1 prop got "
              << fleet.predicted_catch_age_proportion[0][0]
              << " expected " << expected_age1 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.predicted_catch_age_proportion[0][9],
                    expected_age10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: age10 prop got "
              << fleet.predicted_catch_age_proportion[0][9]
              << " expected " << expected_age10 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.predicted_catch_age_proportion[1][9],
                    expected_age10)) {
    std::cerr << "FAIL: scaled year should have same proportions\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA catch age composition observation regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_12_observation_agecomp_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/12_observation_agecomp_caa/bigeye_v2_12_observation_agecomp_caa_check.cpp \
  -o build/examples/bigeye_v2_12_observation_agecomp_caa_check

./build/examples/bigeye_v2_12_observation_agecomp_caa_check
SH

chmod +x run_bigeye_v2_12_observation_agecomp_caa_check.sh

echo "migrated CAA catch age composition observation step"
