#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level04_catch_check"

cat > "$BASE/common/catch.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace bigeye_v2 {

struct BaranovCatch {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> & /*p*/,
                  BigeyeDerived<T> &d) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    d.catch_numbers_at_age.assign(ny, std::array<T, kAges>{});
    d.catch_biomass_at_age.assign(ny, std::array<T, kAges>{});
    d.total_catch_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        const T z = d.z_at_age[a];

        const T catch_numbers =
            d.numbers_at_age[y][a] *
            (d.f_at_age[a] / z) *
            (T(1.0) - std::exp(-z));

        d.catch_numbers_at_age[y][a] = catch_numbers;
        d.catch_biomass_at_age[y][a] = catch_numbers * d.weight_at_age[a];
        d.total_catch_biomass_by_year[y] += d.catch_biomass_at_age[y][a];
      }
    }
  }
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/derived.hpp")
s = p.read_text()
s = s.replace(
"""  std::vector<T> spawning_biomass_by_year{};
};""",
"""  std::vector<T> spawning_biomass_by_year{};

  std::vector<std::array<T, kAges>> catch_numbers_at_age{};
  std::vector<std::array<T, kAges>> catch_biomass_at_age{};
  std::vector<T> total_catch_biomass_by_year{};
};"""
)
p.write_text(s)
PY

cat > "$BASE/level04_catch_check/bigeye_v2_level04_catch_check.cpp" <<'CPP'
#include "../common/catch.hpp"
#include "../common/life_history.hpp"
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

  constexpr double expected_catch_n_age_1 = 3.041727836159226;
  constexpr double expected_catch_n_age_5 = 45.15801288455608;
  constexpr double expected_catch_n_age_10 = 15.49860987509641;
  constexpr double expected_total_catch_biomass = 2113.704431353481;

  for (int y = 0; y < data.n_years; ++y) {
    if (!nearly_equal(d.catch_numbers_at_age[y][0], expected_catch_n_age_1)) {
      std::cerr << std::setprecision(17)
                << "FAIL: catch_numbers_at_age[" << y << "][0] got "
                << d.catch_numbers_at_age[y][0]
                << " expected " << expected_catch_n_age_1
                << " diff "
                << (d.catch_numbers_at_age[y][0] - expected_catch_n_age_1)
                << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.catch_numbers_at_age[y][4], expected_catch_n_age_5)) {
      std::cerr << std::setprecision(17)
                << "FAIL: catch_numbers_at_age[" << y << "][4] got "
                << d.catch_numbers_at_age[y][4]
                << " expected " << expected_catch_n_age_5
                << " diff "
                << (d.catch_numbers_at_age[y][4] - expected_catch_n_age_5)
                << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.catch_numbers_at_age[y][9], expected_catch_n_age_10)) {
      std::cerr << std::setprecision(17)
                << "FAIL: catch_numbers_at_age[" << y << "][9] got "
                << d.catch_numbers_at_age[y][9]
                << " expected " << expected_catch_n_age_10
                << " diff "
                << (d.catch_numbers_at_age[y][9] - expected_catch_n_age_10)
                << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.total_catch_biomass_by_year[y],
                      expected_total_catch_biomass)) {
      std::cerr << std::setprecision(17)
                << "FAIL: total_catch_biomass_by_year[" << y << "] got "
                << d.total_catch_biomass_by_year[y]
                << " expected " << expected_total_catch_biomass
                << " diff "
                << (d.total_catch_biomass_by_year[y] -
                    expected_total_catch_biomass)
                << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 Level04 catch regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level04_catch_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level04_catch_check/bigeye_v2_level04_catch_check.cpp \
  -o build/examples/bigeye_v2_level04_catch_check

./build/examples/bigeye_v2_level04_catch_check
SH

chmod +x run_bigeye_v2_level04_catch_check.sh

echo "created Bigeye v2 Level04 catch check"
