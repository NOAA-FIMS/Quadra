#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level03_mortality_check"

cat > "$BASE/common/mortality.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

namespace bigeye_v2 {

struct FishingMortality {
  template <typename T>
  void operator()(const BigeyeModelData<T> & /*data*/,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    for (int a = 0; a < kAges; ++a) {
      d.f_at_age[a] = p.fbar * d.selectivity_at_age[a];
      d.z_at_age[a] = d.m_at_age[a] + d.f_at_age[a];
    }
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/model_parameters.hpp")
s = p.read_text()
s = s.replace(
"""  T sel_slope = T(1.0);""",
"""  T sel_slope = T(1.0);
  T fbar = T(0.2);"""
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/derived.hpp")
s = p.read_text()
s = s.replace(
"""  std::array<T, kAges> selectivity_at_age{};""",
"""  std::array<T, kAges> selectivity_at_age{};
  std::array<T, kAges> f_at_age{};
  std::array<T, kAges> z_at_age{};"""
)
p.write_text(s)
PY

cat > "$BASE/level03_mortality_check/bigeye_v2_level03_mortality_check.cpp" <<'CPP'
#include "../common/life_history.hpp"
#include "../common/selectivity.hpp"
#include "../common/mortality.hpp"

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

  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.fbar = 0.2;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  LogisticSelectivity{}(data, p, d);
  FishingMortality{}(data, p, d);

  constexpr double expected_f[kAges] = {
      0.003597241992418312,
      0.009485174635513357,
      0.023840584404423512,
      0.053788284273999020,
      0.100000000000000000,
      0.146211715726000980,
      0.176159415595576460,
      0.190514825364486680,
      0.196402758007581700,
      0.198661429815143070};

  constexpr double expected_z[kAges] = {
      0.341097241992418300,
      0.346985174635513300,
      0.361340584404423500,
      0.503788284273999000,
      0.550000000000000000,
      0.596211715726001000,
      0.626159415595576500,
      0.483014825364486700,
      0.488902758007581700,
      0.491161429815143100};

  for (int a = 0; a < kAges; ++a) {
    if (!nearly_equal(d.f_at_age[a], expected_f[a])) {
      std::cerr << std::setprecision(17)
                << "FAIL: f_at_age[" << a << "] got " << d.f_at_age[a]
                << " expected " << expected_f[a]
                << " diff " << (d.f_at_age[a] - expected_f[a]) << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.z_at_age[a], expected_z[a])) {
      std::cerr << std::setprecision(17)
                << "FAIL: z_at_age[" << a << "] got " << d.z_at_age[a]
                << " expected " << expected_z[a]
                << " diff " << (d.z_at_age[a] - expected_z[a]) << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 Level03 mortality regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level03_mortality_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level03_mortality_check/bigeye_v2_level03_mortality_check.cpp \
  -o build/examples/bigeye_v2_level03_mortality_check

./build/examples/bigeye_v2_level03_mortality_check
SH

chmod +x run_bigeye_v2_level03_mortality_check.sh

echo "created Bigeye v2 Level03 mortality check"
