#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

cat > "$BASE/architecture/steps/life_history/life_history.hpp" <<'CPP'
#pragma once

#include "../../../common/bigeye_constants.hpp"
#include "../../parameters/life_history_parameters.hpp"
#include "../../state/life_history_state.hpp"
#include "../../../common/model_data.hpp"

#include <cmath>

namespace bigeye_v2 {

// Bigeye life history.
// Computes natural mortality, weight, and maturity at age.
struct BigeyeLifeHistory {
  template <typename T>
  void operator()(const BigeyeModelData<T> &,
                  const LifeHistoryParameters<T> &p,
                  LifeHistoryState<T> &life) const {
    const T adult_m = T(0.45);
    const T m_young = adult_m * std::exp(p.log_m_young_offset);
    const T m_old = adult_m * std::exp(p.log_m_old_offset);

    for (int a = 0; a < kAges; ++a) {
      const int age = a + 1;

      if (age <= 3) {
        life.m_at_age[a] = m_young;
      } else if (age >= 8) {
        life.m_at_age[a] = m_old;
      } else {
        life.m_at_age[a] = adult_m;
      }

      life.weight_at_age[a] = T(2 * age - 1);
      life.maturity_at_age[a] = age >= 4 ? T(1.0) : T(0.0);
    }
  }
};

} // namespace bigeye_v2
CPP

mkdir -p "$BASE/01_life_history_caa"

cat > "$BASE/01_life_history_caa/bigeye_v2_01_life_history_caa_check.cpp" <<'CPP'
#include "../architecture/steps/life_history/life_history.hpp"

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

  LifeHistoryParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);

  LifeHistoryState<double> life;

  BigeyeLifeHistory{}(data, p, life);

  constexpr double expected_m_age_1 = 0.3375;
  constexpr double expected_m_age_4 = 0.45;
  constexpr double expected_m_age_8 = 0.2925;

  if (!nearly_equal(life.m_at_age[0], expected_m_age_1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: m age 1 got " << life.m_at_age[0] << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.m_at_age[3], expected_m_age_4)) {
    std::cerr << std::setprecision(17)
              << "FAIL: m age 4 got " << life.m_at_age[3] << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.m_at_age[7], expected_m_age_8)) {
    std::cerr << std::setprecision(17)
              << "FAIL: m age 8 got " << life.m_at_age[7] << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.weight_at_age[9], 19.0)) {
    std::cerr << "FAIL: weight age 10\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.maturity_at_age[2], 0.0) ||
      !nearly_equal(life.maturity_at_age[3], 1.0)) {
    std::cerr << "FAIL: maturity knife-edge\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA life-history regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_01_life_history_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/01_life_history_caa/bigeye_v2_01_life_history_caa_check.cpp \
  -o build/examples/bigeye_v2_01_life_history_caa_check

./build/examples/bigeye_v2_01_life_history_caa_check
SH

chmod +x run_bigeye_v2_01_life_history_caa_check.sh

echo "migrated CAA life-history step"
