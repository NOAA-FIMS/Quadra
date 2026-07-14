#include "../architecture/packages/observation/observation_package.hpp"

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

  BigeyeModelData<double> data;
  data.n_years = 2;

  FleetParameters<double> parameters;
  parameters.q_index = 0.01;

  PopulationState<double> population;
  population.spawning_biomass_by_year = {10000.0, 12000.0};

  FleetState<double> fleet;
  fleet.catch_numbers_at_age.assign(static_cast<std::size_t>(data.n_years),
                                    std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    fleet.catch_numbers_at_age[0][a] = static_cast<double>(a + 1);
    fleet.catch_numbers_at_age[1][a] = 2.0 * static_cast<double>(a + 1);
  }

  ObservationPackage{}(data, parameters, population, fleet);

  if (!nearly_equal(fleet.predicted_index_by_year[0], 100.0) ||
      !nearly_equal(fleet.predicted_index_by_year[1], 120.0)) {
    std::cerr << "FAIL: biomass index prediction\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.predicted_catch_age_proportion[0][0], 1.0 / 55.0) ||
      !nearly_equal(fleet.predicted_catch_age_proportion[0][9], 10.0 / 55.0)) {
    std::cerr << std::setprecision(17)
              << "FAIL: age composition prediction got age1 "
              << fleet.predicted_catch_age_proportion[0][0] << " age10 "
              << fleet.predicted_catch_age_proportion[0][9] << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA ObservationPackage regression\n";
  return EXIT_SUCCESS;
}
