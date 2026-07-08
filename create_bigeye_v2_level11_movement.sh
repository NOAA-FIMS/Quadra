#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/level11_movement_check"

cat > "$BASE/common/movement.hpp" <<'CPP'
#pragma once

#include "bigeye_constants.hpp"
#include "model_data.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T>
struct BigeyePopulationDerived {
  std::vector<std::array<T, kAges>> numbers_at_age{};
};

template <typename T>
struct BigeyeMovementParameters {
  T move_0_to_1 = T(0.10);
  T move_1_to_0 = T(0.05);
};

struct TwoPopulationMovement {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeMovementParameters<T> &p,
                  std::vector<BigeyePopulationDerived<T>> &pop) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    if (pop.size() != 2) {
      return;
    }

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        const T n0 = pop[0].numbers_at_age[y][a];
        const T n1 = pop[1].numbers_at_age[y][a];

        pop[0].numbers_at_age[y][a] =
            n0 * (T(1.0) - p.move_0_to_1) + n1 * p.move_1_to_0;

        pop[1].numbers_at_age[y][a] =
            n1 * (T(1.0) - p.move_1_to_0) + n0 * p.move_0_to_1;
      }
    }
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/level11_movement_check/bigeye_v2_level11_movement_check.cpp" <<'CPP'
#include "../common/movement.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

bool nearly_equal(double a, double b, double tol = 1.0e-9) {
  return std::abs(a - b) <= tol;
}

} // namespace

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 2;

  std::vector<BigeyePopulationDerived<double>> pop(2);

  for (auto &x : pop) {
    x.numbers_at_age.assign(
        static_cast<std::size_t>(data.n_years),
        std::array<double, kAges>{});
  }

  // Population 0 starts larger.
  // Population 1 starts smaller.
  for (int a = 0; a < kAges; ++a) {
    pop[0].numbers_at_age[0][a] = 1000.0 + 10.0 * a;
    pop[1].numbers_at_age[0][a] = 500.0 + 5.0 * a;

    pop[0].numbers_at_age[1][a] = 2000.0 + 20.0 * a;
    pop[1].numbers_at_age[1][a] = 1000.0 + 10.0 * a;
  }

  BigeyeMovementParameters<double> mp;
  mp.move_0_to_1 = 0.10;
  mp.move_1_to_0 = 0.05;

  TwoPopulationMovement{}(data, mp, pop);

  constexpr double expected_pop0_year0_age1 = 925.0;
  constexpr double expected_pop1_year0_age1 = 575.0;

  constexpr double expected_pop0_year1_age10 = 1932.5;
  constexpr double expected_pop1_year1_age10 = 1202.5;

  if (!nearly_equal(pop[0].numbers_at_age[0][0],
                    expected_pop0_year0_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: pop0 y0 age1 got "
              << pop[0].numbers_at_age[0][0]
              << " expected " << expected_pop0_year0_age1 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop[1].numbers_at_age[0][0],
                    expected_pop1_year0_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: pop1 y0 age1 got "
              << pop[1].numbers_at_age[0][0]
              << " expected " << expected_pop1_year0_age1 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop[0].numbers_at_age[1][9],
                    expected_pop0_year1_age10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: pop0 y1 age10 got "
              << pop[0].numbers_at_age[1][9]
              << " expected " << expected_pop0_year1_age10 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop[1].numbers_at_age[1][9],
                    expected_pop1_year1_age10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: pop1 y1 age10 got "
              << pop[1].numbers_at_age[1][9]
              << " expected " << expected_pop1_year1_age10 << "\n";
    return EXIT_FAILURE;
  }

  const double before_total_y0_age1 = 1000.0 + 500.0;
  const double after_total_y0_age1 =
      pop[0].numbers_at_age[0][0] + pop[1].numbers_at_age[0][0];

  if (!nearly_equal(after_total_y0_age1, before_total_y0_age1)) {
    std::cerr << "FAIL: movement should conserve total abundance\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level11 movement regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_level11_movement_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level11_movement_check/bigeye_v2_level11_movement_check.cpp \
  -o build/examples/bigeye_v2_level11_movement_check

./build/examples/bigeye_v2_level11_movement_check
SH

chmod +x run_bigeye_v2_level11_movement_check.sh

echo "created Bigeye v2 Level11 movement check"
