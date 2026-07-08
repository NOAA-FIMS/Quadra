#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level05_catch_likelihood_check"

cat > "$BASE/common/likelihood.hpp" <<'CPP'
#pragma once

#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LognormalCatchLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.catch_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_catch_biomass_by_year.size(); ++y) {
      const T obs = data.observed_catch_biomass_by_year[y];
      const T pred = d.total_catch_biomass_by_year[y];

      const T r = (std::log(obs) - std::log(pred)) / p.catch_sigma;

      d.catch_nll +=
          T(0.5) * r * r +
          std::log(p.catch_sigma) +
          std::log(obs);
    }

    d.total_nll += d.catch_nll;
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/model_data.hpp")
s = p.read_text()
s = s.replace(
"""#pragma once

namespace bigeye_v2 {""",
"""#pragma once

#include <vector>

namespace bigeye_v2 {"""
)
s = s.replace(
"""  int n_years = 1;""",
"""  int n_years = 1;
  std::vector<T> observed_catch_biomass_by_year{};"""
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/model_parameters.hpp")
s = p.read_text()
s = s.replace(
"""  T fbar = T(0.2);""",
"""  T fbar = T(0.2);
  T catch_sigma = T(0.1);"""
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/derived.hpp")
s = p.read_text()
s = s.replace(
"""  std::vector<T> total_catch_biomass_by_year{};
};""",
"""  std::vector<T> total_catch_biomass_by_year{};

  T catch_nll = T(0.0);
  T total_nll = T(0.0);
};"""
)
p.write_text(s)
PY

cat > "$BASE/level05_catch_likelihood_check/bigeye_v2_level05_catch_likelihood_check.cpp" <<'CPP'
#include "../common/catch.hpp"
#include "../common/life_history.hpp"
#include "../common/likelihood.hpp"
#include "../common/mortality.hpp"
#include "../common/population.hpp"
#include "../common/recruitment.hpp"
#include "../common/selectivity.hpp"

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
  data.observed_catch_biomass_by_year = {
      1326.1639007786976,
      1350.0,
      1300.0};

  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);
  p.r0 = 1000.0;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.fbar = 0.2;
  p.catch_sigma = 0.1;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  LogisticSelectivity{}(data, p, d);
  FishingMortality{}(data, p, d);
  UnfishedPopulation{}(data, p, d);
  BaranovCatch{}(data, p, d);
  LognormalCatchLikelihood{}(data, p, d);

  constexpr double expected_catch_nll = 19.18886774409258;
  constexpr double expected_total_nll = 19.18886774409258;

  if (!nearly_equal(d.catch_nll, expected_catch_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch_nll got " << d.catch_nll
              << " expected " << expected_catch_nll
              << " diff " << (d.catch_nll - expected_catch_nll) << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.total_nll, expected_total_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: total_nll got " << d.total_nll
              << " expected " << expected_total_nll
              << " diff " << (d.total_nll - expected_total_nll) << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level05 catch likelihood regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level05_catch_likelihood_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level05_catch_likelihood_check/bigeye_v2_level05_catch_likelihood_check.cpp \
  -o build/examples/bigeye_v2_level05_catch_likelihood_check

./build/examples/bigeye_v2_level05_catch_likelihood_check
SH

chmod +x run_bigeye_v2_level05_catch_likelihood_check.sh

echo "created Bigeye v2 Level05 catch likelihood check"
