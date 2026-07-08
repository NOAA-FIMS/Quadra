#include "../architecture/state/fleet_state.hpp"
#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/population/survival.hpp"

#include <array>
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

  PopulationState<double> population;
  population.numbers_at_age.assign(2, std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    population.numbers_at_age[0][a] = 1000.0;
    population.numbers_at_age[1][a] = 2000.0;
  }

  FleetState<double> fleet;
  for (int a = 0; a < kAges; ++a) {
    fleet.z_at_age[a] = 0.1 * static_cast<double>(a + 1);
  }

  Survival{}(fleet, population);

  constexpr double expected_y0_age1 = 904.83741803595957;
  constexpr double expected_y0_age10 = 367.87944117144235;
  constexpr double expected_y1_age1 = 1809.6748360719191;

  if (!nearly_equal(population.survivors_at_age[0][0], expected_y0_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y0 age1 survivor got "
              << population.survivors_at_age[0][0]
              << " expected " << expected_y0_age1
              << " diff " << (population.survivors_at_age[0][0] -
                               expected_y0_age1)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.survivors_at_age[0][9], expected_y0_age10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y0 age10 survivor got "
              << population.survivors_at_age[0][9]
              << " expected " << expected_y0_age10
              << " diff " << (population.survivors_at_age[0][9] -
                               expected_y0_age10)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(population.survivors_at_age[1][0], expected_y1_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age1 survivor got "
              << population.survivors_at_age[1][0]
              << " expected " << expected_y1_age1
              << " diff " << (population.survivors_at_age[1][0] -
                               expected_y1_age1)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA population survival regression\n";
  return EXIT_SUCCESS;
}
