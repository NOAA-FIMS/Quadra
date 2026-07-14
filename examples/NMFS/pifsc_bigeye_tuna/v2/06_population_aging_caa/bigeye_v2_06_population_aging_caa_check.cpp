#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/population/aging.hpp"

#include <array>
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

  PopulationState<double> population;
  population.survivors_at_age.assign(3, std::array<double, kAges>{});
  population.numbers_at_age.assign(3, std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    population.survivors_at_age[0][a] = 100.0 + a;
    population.survivors_at_age[1][a] = 200.0 + a;
  }

  Aging{}(population);

  constexpr double expected_y1_age2 = 100.0;
  constexpr double expected_y1_age9 = 107.0;
  constexpr double expected_y2_age2 = 200.0;

  if (!nearly_equal(population.numbers_at_age[1][1], expected_y1_age2)) {
    std::cerr << std::setprecision(17) << "FAIL: y1 age2 got "
              << population.numbers_at_age[1][1] << " expected "
              << expected_y1_age2 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.numbers_at_age[1][8], expected_y1_age9)) {
    std::cerr << std::setprecision(17) << "FAIL: y1 age9 got "
              << population.numbers_at_age[1][8] << " expected "
              << expected_y1_age9 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.numbers_at_age[2][1], expected_y2_age2)) {
    std::cerr << std::setprecision(17) << "FAIL: y2 age2 got "
              << population.numbers_at_age[2][1] << " expected "
              << expected_y2_age2 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.numbers_at_age[1][0], 0.0)) {
    std::cerr << "FAIL: aging should not set recruits/age1\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.numbers_at_age[1][9], 0.0)) {
    std::cerr << "FAIL: aging should not set plus group\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA population aging regression\n";
  return EXIT_SUCCESS;
}
