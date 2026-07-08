#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

cat > "$BASE/architecture/steps/fleet/logistic_selectivity.hpp" <<'CPP'
#pragma once

#include "../../../common/bigeye_constants.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../../common/model_data.hpp"

#include <cmath>

namespace bigeye_v2 {

// Logistic selectivity.
// Computes vulnerability-at-age for one fleet.
struct LogisticSelectivity {
  template <typename T>
  void operator()(const BigeyeModelData<T> &,
                  const FleetParameters<T> &p,
                  FleetState<T> &fleet) const {
    for (int a = 0; a < kAges; ++a) {
      const T age = T(a + 1);
      fleet.selectivity_at_age[a] =
          T(1.0) / (T(1.0) + std::exp(-p.sel_slope * (age - p.sel_a50)));
    }
  }
};

} // namespace bigeye_v2
CPP

mkdir -p "$BASE/level02_selectivity_caa_check"

cat > "$BASE/level02_selectivity_caa_check/bigeye_v2_level02_selectivity_caa_check.cpp" <<'CPP'
#include "../architecture/steps/fleet/logistic_selectivity.hpp"

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

  BigeyeModelData<double> data;

  FleetParameters<double> p;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;

  FleetState<double> fleet;

  LogisticSelectivity{}(data, p, fleet);

  constexpr double expected[kAges] = {
      0.01798620996209156,
      0.04742587317756678,
      0.11920292202211755,
      0.26894142136999510,
      0.50000000000000000,
      0.73105857863000490,
      0.88079707797788230,
      0.95257412682243340,
      0.98201379003790850,
      0.99330714907571530};

  for (int a = 0; a < kAges; ++a) {
    if (!nearly_equal(fleet.selectivity_at_age[a], expected[a])) {
      std::cerr << std::setprecision(17)
                << "FAIL: CAA selectivity_at_age[" << a << "] got "
                << fleet.selectivity_at_age[a]
                << " expected " << expected[a]
                << " diff " << (fleet.selectivity_at_age[a] - expected[a])
                << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 CAA logistic selectivity regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level02_selectivity_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level02_selectivity_caa_check/bigeye_v2_level02_selectivity_caa_check.cpp \
  -o build/examples/bigeye_v2_level02_selectivity_caa_check

./build/examples/bigeye_v2_level02_selectivity_caa_check
SH

chmod +x run_bigeye_v2_level02_selectivity_caa_check.sh

echo "migrated CAA logistic selectivity step"
