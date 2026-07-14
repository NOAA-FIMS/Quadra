#include "../architecture/state/fleet_state.hpp"
#include "../architecture/steps/observation/catch_age_composition.hpp"
#include "../common/model_data.hpp"

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

  FleetState<double> fleet;
  fleet.catch_numbers_at_age.assign(static_cast<std::size_t>(data.n_years),
                                    std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    fleet.catch_numbers_at_age[0][a] = static_cast<double>(a + 1);
    fleet.catch_numbers_at_age[1][a] = 2.0 * static_cast<double>(a + 1);
  }

  CatchAgeCompositionPrediction{}(data, fleet);

  constexpr double expected_age1 = 1.0 / 55.0;
  constexpr double expected_age10 = 10.0 / 55.0;

  if (!nearly_equal(fleet.predicted_catch_age_proportion[0][0],
                    expected_age1)) {
    std::cerr << std::setprecision(17) << "FAIL: age1 prop got "
              << fleet.predicted_catch_age_proportion[0][0] << " expected "
              << expected_age1 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.predicted_catch_age_proportion[0][9],
                    expected_age10)) {
    std::cerr << std::setprecision(17) << "FAIL: age10 prop got "
              << fleet.predicted_catch_age_proportion[0][9] << " expected "
              << expected_age10 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.predicted_catch_age_proportion[1][9],
                    expected_age10)) {
    std::cerr << "FAIL: scaled year should have same proportions\n";
    return EXIT_FAILURE;
  }

  std::cout
      << "PASSED: Bigeye v2 CAA catch age composition observation regression\n";
  return EXIT_SUCCESS;
}
