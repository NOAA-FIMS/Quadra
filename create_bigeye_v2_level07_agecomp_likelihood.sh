#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level07_agecomp_likelihood_check"

cat > "$BASE/common/agecomp.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

namespace bigeye_v2 {

struct CatchAgeCompositionPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &,
                  BigeyeDerived<T> &d) const {
    const auto ny = static_cast<std::size_t>(data.n_years);
    d.predicted_catch_age_proportion.assign(ny, std::array<T, kAges>{});

    for (std::size_t y = 0; y < ny; ++y) {
      T total = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        total += d.catch_numbers_at_age[y][a];
      }

      for (int a = 0; a < kAges; ++a) {
        d.predicted_catch_age_proportion[y][a] =
            d.catch_numbers_at_age[y][a] / total;
      }
    }
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/model_data.hpp")
s = p.read_text()
s = s.replace("#include <vector>", "#include \"bigeye_constants.hpp\"\n\n#include <array>\n#include <vector>")
s = s.replace(
"  std::vector<T> observed_index_by_year{};",
"  std::vector<T> observed_index_by_year{};\n  std::vector<std::array<T, kAges>> observed_catch_age_proportion{};\n  std::vector<T> catch_agecomp_sample_size{};"
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/derived.hpp")
s = p.read_text()
s = s.replace(
"  std::vector<T> predicted_index_by_year{};",
"  std::vector<T> predicted_index_by_year{};\n  std::vector<std::array<T, kAges>> predicted_catch_age_proportion{};"
)
s = s.replace(
"  T index_nll = T(0.0);",
"  T index_nll = T(0.0);\n  T agecomp_nll = T(0.0);"
)
p.write_text(s)
PY

python3 - <<'PY'
from pathlib import Path
p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/likelihood.hpp")
s = p.read_text()

insert = r'''

struct MultinomialAgeCompLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &,
                  BigeyeDerived<T> &d) const {
    d.agecomp_nll = T(0.0);

    const T eps = T(1.0e-12);

    for (std::size_t y = 0; y < data.observed_catch_age_proportion.size(); ++y) {
      const T n_eff = data.catch_agecomp_sample_size[y];

      for (int a = 0; a < kAges; ++a) {
        const T obs = data.observed_catch_age_proportion[y][a];
        const T pred = d.predicted_catch_age_proportion[y][a] + eps;

        d.agecomp_nll -= n_eff * obs * std::log(pred);
      }
    }

    d.total_nll += d.agecomp_nll;
  }
};
'''

s = s.replace("\n} // namespace bigeye_v2\n", insert + "\n} // namespace bigeye_v2\n")
p.write_text(s)
PY

cat > "$BASE/level07_agecomp_likelihood_check/bigeye_v2_level07_agecomp_likelihood_check.cpp" <<'CPP'
#include "../common/agecomp.hpp"
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
  data.catch_agecomp_sample_size = {100.0, 100.0, 100.0};

  data.observed_catch_age_proportion = {
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06}};

  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);
  p.r0 = 1000.0;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.fbar = 0.2;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  LogisticSelectivity{}(data, p, d);
  FishingMortality{}(data, p, d);
  UnfishedPopulation{}(data, p, d);
  BaranovCatch{}(data, p, d);
  CatchAgeCompositionPrediction{}(data, p, d);
  MultinomialAgeCompLikelihood{}(data, p, d);

  constexpr double expected_prop_age_1 = 0.016592518658389734;
  constexpr double expected_prop_age_10 = 0.11272888052337992;
  constexpr double expected_agecomp_nll = 684.1893174267169;

  if (!nearly_equal(d.predicted_catch_age_proportion[0][0], expected_prop_age_1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: prop age 1 got "
              << d.predicted_catch_age_proportion[0][0]
              << " expected " << expected_prop_age_1
              << " diff "
              << (d.predicted_catch_age_proportion[0][0] - expected_prop_age_1)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.predicted_catch_age_proportion[0][9], expected_prop_age_10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: prop age 10 got "
              << d.predicted_catch_age_proportion[0][9]
              << " expected " << expected_prop_age_10
              << " diff "
              << (d.predicted_catch_age_proportion[0][9] - expected_prop_age_10)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.agecomp_nll, expected_agecomp_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: agecomp_nll got " << d.agecomp_nll
              << " expected " << expected_agecomp_nll
              << " diff " << (d.agecomp_nll - expected_agecomp_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level07 age composition likelihood regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level07_agecomp_likelihood_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level07_agecomp_likelihood_check/bigeye_v2_level07_agecomp_likelihood_check.cpp \
  -o build/examples/bigeye_v2_level07_agecomp_likelihood_check

./build/examples/bigeye_v2_level07_agecomp_likelihood_check
SH

chmod +x run_bigeye_v2_level07_agecomp_likelihood_check.sh

echo "created Bigeye v2 Level07 age composition likelihood check"
