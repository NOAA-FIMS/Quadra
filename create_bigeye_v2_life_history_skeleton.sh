#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/common" "$BASE/level00_life_history_check"

cat > "$BASE/common/bigeye_constants.hpp" <<'CPP'
#pragma once

namespace bigeye_v2 {

inline constexpr int kAges = 10;

} // namespace bigeye_v2
CPP

cat > "$BASE/common/life_history.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"

#include <array>
#include <cmath>

namespace bigeye_v2 {

template <typename T>
struct BigeyeLifeHistory {
  // Parameters
  T log_m_young_offset = T(0.0);
  T log_m_old_offset = T(0.0);

  // Derived quantities
  std::array<T, kAges> m_at_age{};
  std::array<T, kAges> weight_at_age{};
  std::array<T, kAges> maturity_at_age{};

  void operator()() {
    const T adult_m = T(0.45);
    const T m_young = adult_m * std::exp(log_m_young_offset);
    const T m_old = adult_m * std::exp(log_m_old_offset);

    for (int a = 0; a < kAges; ++a) {
      if (a <= 2) {
        m_at_age[a] = m_young;
      } else if (a >= 7) {
        m_at_age[a] = m_old;
      } else {
        m_at_age[a] = adult_m;
      }

      // Placeholder deterministic curves.
      // Replace later with Bigeye-specific values.
      weight_at_age[a] = T(1.0 + 2.0 * a);
      maturity_at_age[a] = T(a >= 3 ? 1.0 : 0.0);
    }
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/level00_life_history_check/bigeye_v2_level00_life_history_check.cpp" <<'CPP'
#include "../common/life_history.hpp"

#include <iomanip>
#include <iostream>

int main() {
  using namespace bigeye_v2;

  BigeyeLifeHistory<double> lh;
  lh.log_m_young_offset = std::log(0.75);
  lh.log_m_old_offset = std::log(0.65);

  lh();

  std::cout << "age,m_at_age,weight_at_age,maturity_at_age\n";
  for (int a = 0; a < kAges; ++a) {
    std::cout << (a + 1) << ","
              << std::setprecision(12) << lh.m_at_age[a] << ","
              << lh.weight_at_age[a] << ","
              << lh.maturity_at_age[a] << "\n";
  }

  return 0;
}
CPP

cat > run_bigeye_v2_level00_life_history_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level00_life_history_check/bigeye_v2_level00_life_history_check.cpp \
  -o build/examples/bigeye_v2_level00_life_history_check

./build/examples/bigeye_v2_level00_life_history_check
SH

chmod +x run_bigeye_v2_level00_life_history_check.sh

echo "Created Bigeye v2 life-history skeleton."
echo "Run:"
echo "  ./run_bigeye_v2_level00_life_history_check.sh"
