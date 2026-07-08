#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level06_index_likelihood_check"

cat > "$BASE/common/index.hpp" <<'CPP'
#pragma once

#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

namespace bigeye_v2 {

struct BiomassIndexPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.predicted_index_by_year.assign(
        static_cast<std::size_t>(data.n_years), T(0.0));

    for (std::size_t y = 0; y < d.spawning_biomass_by_year.size(); ++y) {
      d.predicted_index_by_year[y] = p.q_index * d.spawning_biomass_by_year[y];
    }
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/model_data.hpp")
s = p.read_text()
s = s.replace(
"  std::vector<T> observed_catch_biomass_by_year{};",
"  std::vector<T> observed_catch_biomass_by_year{};\n  std::vector<T> observed_index_by_year{};"
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/model_parameters.hpp")
s = p.read_text()
s = s.replace(
"  T catch_sigma = T(0.1);",
"  T catch_sigma = T(0.1);\n  T q_index = T(0.01);\n  T index_sigma = T(0.2);"
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/derived.hpp")
s = p.read_text()
s = s.replace(
"  std::vector<T> total_catch_biomass_by_year{};",
"  std::vector<T> total_catch_biomass_by_year{};\n  std::vector<T> predicted_index_by_year{};"
)
s = s.replace(
"  T catch_nll = T(0.0);",
"  T catch_nll = T(0.0);\n  T index_nll = T(0.0);"
)
p.write_text(s)
PY

python3 - <<'PY'
from pathlib import Path
p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/likelihood.hpp")
s = p.read_text()
insert = r'''

struct LognormalIndexLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.index_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_index_by_year.size(); ++y) {
      const T obs = data.observed_index_by_year[y];
      const T pred = d.predicted_index_by_year[y];

      const T r = (std::log(obs) - std::log(pred)) / p.index_sigma;

      d.index_nll +=
          T(0.5) * r * r +
          std::log(p.index_sigma) +
          std::log(obs);
    }

    d.total_nll += d.index_nll;
  }
};
'''
s = s.replace("\n} // namespace bigeye_v2\n", insert + "\n} // namespace bigeye_v2\n")
p.write_text(s)
PY

cat > "$BASE/level06_index_likelihood_check/bigeye_v2_level06_index_likelihood_check.cpp" <<'CPP'
#include "../common/catch.hpp"
#include "../common/index.hpp"
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
  data.observed_index_by_year = {116.46701723019194, 120.0, 110.0};

  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);
  p.r0 = 1000.0;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.fbar = 0.2;
  p.q_index = 0.01;
  p.index_sigma = 0.2;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  LogisticSelectivity{}(data, p, d);
  FishingMortality{}(data, p, d);
  UnfishedPopulation{}(data, p, d);
  BiomassIndexPrediction{}(data, p, d);
  LognormalIndexLikelihood{}(data, p, d);

  constexpr double expected_pred_index = 116.46701723019194;
  constexpr double expected_index_nll = 9.858578377283694;

  if (!nearly_equal(d.predicted_index_by_year[0], expected_pred_index)) {
    std::cerr << std::setprecision(17)
              << "FAIL: predicted_index got " << d.predicted_index_by_year[0]
              << " expected " << expected_pred_index << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.index_nll, expected_index_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: index_nll got " << d.index_nll
              << " expected " << expected_index_nll
              << " diff " << (d.index_nll - expected_index_nll) << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level06 index likelihood regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level06_index_likelihood_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level06_index_likelihood_check/bigeye_v2_level06_index_likelihood_check.cpp \
  -o build/examples/bigeye_v2_level06_index_likelihood_check

./build/examples/bigeye_v2_level06_index_likelihood_check
SH

chmod +x run_bigeye_v2_level06_index_likelihood_check.sh

echo "created Bigeye v2 Level06 index likelihood check"
