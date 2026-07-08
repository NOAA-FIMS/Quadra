#include "../architecture/parameters/fleet_parameters.hpp"
#include "../architecture/state/fleet_state.hpp"
#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/observation/biomass_index.hpp"
#include "../common/model_data.hpp"

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
  data.n_years = 3;

  PopulationState<double> population;
  population.spawning_biomass_by_year = {10000.0, 12000.0, 9000.0};

  FleetParameters<double> parameters;
  parameters.q_index = 0.01;

  FleetState<double> fleet;

  BiomassIndexPrediction{}(data, parameters, population, fleet);

  constexpr double expected_y0 = 100.0;
  constexpr double expected_y1 = 120.0;
  constexpr double expected_y2 = 90.0;

  if (!nearly_equal(fleet.predicted_index_by_year[0], expected_y0) ||
      !nearly_equal(fleet.predicted_index_by_year[1], expected_y1) ||
      !nearly_equal(fleet.predicted_index_by_year[2], expected_y2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: predicted index values got "
              << fleet.predicted_index_by_year[0] << ", "
              << fleet.predicted_index_by_year[1] << ", "
              << fleet.predicted_index_by_year[2] << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA biomass index observation regression\n";
  return EXIT_SUCCESS;
}
