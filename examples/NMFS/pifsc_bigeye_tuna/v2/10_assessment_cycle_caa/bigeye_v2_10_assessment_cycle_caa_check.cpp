#include "../architecture/assessment/assessment_cycle.hpp"

#include <array>
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

  AssessmentParameters<double> parameters;
  parameters.life.log_m_young_offset = std::log(0.75);
  parameters.life.log_m_old_offset = std::log(0.65);

  parameters.populations.resize(1);
  parameters.populations[0].r0 = 1000.0;

  parameters.fleets.resize(1);
  parameters.fleets[0].sel_a50 = 5.0;
  parameters.fleets[0].sel_slope = 1.0;
  parameters.fleets[0].fbar = 0.2;

  AssessmentState<double> state;
  state.populations.resize(1);
  state.fleets.resize(1);

  state.populations[0].numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    state.populations[0].numbers_at_age[0][a] = 1000.0;
  }

  AssessmentCycle{}(data, parameters, state);

  constexpr double expected_m_age1 = 0.3375;
  constexpr double expected_pop_y1_age2 = 710.98976678322833;
  constexpr double expected_fleet_z_age5 = 0.55;
  constexpr double expected_total_catch_y0 = 11918.418702057177;

  if (!nearly_equal(state.life.m_at_age[0], expected_m_age1)) {
    std::cerr << "FAIL: life history not computed\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.populations[0].numbers_at_age[1][1],
                    expected_pop_y1_age2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: population y1 age2 got "
              << state.populations[0].numbers_at_age[1][1]
              << " expected " << expected_pop_y1_age2
              << " diff "
              << (state.populations[0].numbers_at_age[1][1] -
                  expected_pop_y1_age2)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.fleets[0].z_at_age[4], expected_fleet_z_age5)) {
    std::cerr << "FAIL: fleet mortality not computed\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.fleets[0].total_catch_biomass_by_year[0],
                    expected_total_catch_y0)) {
    std::cerr << std::setprecision(17)
              << "FAIL: total catch y0 got "
              << state.fleets[0].total_catch_biomass_by_year[0]
              << " expected " << expected_total_catch_y0
              << " diff "
              << (state.fleets[0].total_catch_biomass_by_year[0] -
                  expected_total_catch_y0)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA AssessmentCycle regression\n";
  return EXIT_SUCCESS;
}
