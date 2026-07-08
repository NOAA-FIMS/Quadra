#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level02_selectivity_check"

cat > "$BASE/common/selectivity.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LogisticSelectivity {
  template <typename T>
  void operator()(const BigeyeModelData<T> & /*data*/,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    for (int a = 0; a < kAges; ++a) {
      const T age = T(a + 1);
      d.selectivity_at_age[a] =
          T(1.0) / (T(1.0) + std::exp(-p.sel_slope * (age - p.sel_a50)));
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
"""  T r0 = T(1000.0);""",
"""  T r0 = T(1000.0);
  T sel_a50 = T(5.0);
  T sel_slope = T(1.0);"""
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/derived.hpp")
s = p.read_text()
s = s.replace(
"""  std::array<T, kAges> maturity_at_age{};""",
"""  std::array<T, kAges> maturity_at_age{};
  std::array<T, kAges> selectivity_at_age{};"""
)
p.write_text(s)
PY

cat > "$BASE/level02_selectivity_check/bigeye_v2_level02_selectivity_check.cpp" <<'CPP'
#include "../common/life_history.hpp"
#include "../common/selectivity.hpp"

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
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;

  BigeyeDerived<double> d;

  LogisticSelectivity{}(data, p, d);

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
    if (!nearly_equal(d.selectivity_at_age[a], expected[a])) {
      std::cerr << std::setprecision(17)
                << "FAIL: selectivity_at_age[" << a << "] got "
                << d.selectivity_at_age[a]
                << " expected " << expected[a]
                << " diff " << (d.selectivity_at_age[a] - expected[a])
                << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 Level02 selectivity regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level02_selectivity_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level02_selectivity_check/bigeye_v2_level02_selectivity_check.cpp \
  -o build/examples/bigeye_v2_level02_selectivity_check

./build/examples/bigeye_v2_level02_selectivity_check
SH

chmod +x run_bigeye_v2_level02_selectivity_check.sh

echo "created Bigeye v2 Level02 selectivity check"
