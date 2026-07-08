#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level01_population_check"

cat > "$BASE/common/recruitment.hpp" <<'CPP'
#pragma once

#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <vector>

namespace bigeye_v2 {

struct FixedRecruitment {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.recruits_by_year.assign(static_cast<std::size_t>(data.n_years), p.r0);
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/common/population.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <cmath>
#include <vector>

namespace bigeye_v2 {

struct UnfishedPopulation {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> & /*p*/,
                  BigeyeDerived<T> &d) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    d.numbers_at_age.assign(ny, std::array<T, kAges>{});
    d.spawning_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      d.numbers_at_age[y][0] = d.recruits_by_year[y];

      for (int a = 1; a < kAges - 1; ++a) {
        d.numbers_at_age[y][a] =
            d.numbers_at_age[y][a - 1] * std::exp(-d.m_at_age[a - 1]);
      }

      const int plus = kAges - 1;
      d.numbers_at_age[y][plus] =
          d.numbers_at_age[y][plus - 1] *
          std::exp(-d.m_at_age[plus - 1]) /
          (T(1.0) - std::exp(-d.m_at_age[plus]));

      for (int a = 0; a < kAges; ++a) {
        d.spawning_biomass_by_year[y] +=
            d.numbers_at_age[y][a] *
            d.weight_at_age[a] *
            d.maturity_at_age[a];
      }
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
"""  T log_m_old_offset = T(0.0);""",
"""  T log_m_old_offset = T(0.0);
  T r0 = T(1000.0);"""
)
p.write_text(s)

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/common/derived.hpp")
s = p.read_text()
s = s.replace("#include <array>", "#include <array>\n#include <vector>")
s = s.replace(
"""  std::array<T, kAges> maturity_at_age{};
};""",
"""  std::array<T, kAges> maturity_at_age{};

  std::vector<T> recruits_by_year{};
  std::vector<std::array<T, kAges>> numbers_at_age{};
  std::vector<T> spawning_biomass_by_year{};
};"""
)
p.write_text(s)
PY

cat > "$BASE/level01_population_check/bigeye_v2_level01_population_check.cpp" <<'CPP'
#include "../common/life_history.hpp"
#include "../common/population.hpp"
#include "../common/recruitment.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

bool nearly_equal(double a, double b, double tol = 1.0e-9) {
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

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  UnfishedPopulation{}(data, p, d);

  constexpr double expected_recruit = 1000.0;
  constexpr double expected_n_age_1 = 1000.0;
  constexpr double expected_n_age_2 = 713.195287287;
  constexpr double expected_n_age_10 = 915.641069432;
  constexpr double expected_ssb = 44131.2734498;

  for (int y = 0; y < data.n_years; ++y) {
    if (!nearly_equal(d.recruits_by_year[y], expected_recruit)) {
      std::cerr << "FAIL: recruits_by_year[" << y << "]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.numbers_at_age[y][0], expected_n_age_1)) {
      std::cerr << "FAIL: numbers_at_age[" << y << "][0]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.numbers_at_age[y][1], expected_n_age_2)) {
      std::cerr << "FAIL: numbers_at_age[" << y << "][1]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.numbers_at_age[y][9], expected_n_age_10)) {
      std::cerr << "FAIL: numbers_at_age[" << y << "][9]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.spawning_biomass_by_year[y], expected_ssb)) {
      std::cerr << "FAIL: spawning_biomass_by_year[" << y << "] got "
                << d.spawning_biomass_by_year[y] << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 Level01 population regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level01_population_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level01_population_check/bigeye_v2_level01_population_check.cpp \
  -o build/examples/bigeye_v2_level01_population_check

./build/examples/bigeye_v2_level01_population_check
SH

chmod +x run_bigeye_v2_level01_population_check.sh

echo "created Bigeye v2 Level01 population check"
