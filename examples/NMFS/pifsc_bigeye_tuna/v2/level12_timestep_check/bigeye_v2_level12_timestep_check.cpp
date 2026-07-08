#include "../common/fleet.hpp"
#include "../common/life_history.hpp"
#include "../common/timestep.hpp"

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

  BigeyeDerived<double> life;
  BigeyeLifeHistory{}(data, p, life);

  BigeyeFleetParameters<double> fp;
  fp.sel_a50 = 5.0;
  fp.sel_slope = 1.0;
  fp.fbar = 0.2;

  BigeyeFleetDerived<double> fleet;
  FleetLogisticSelectivity{}(data, fp, fleet);
  FleetFishingMortality{}(data, fp, life, fleet);

  TimeStepPopulationState<double> pop;
  pop.numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  // Initial year numbers-at-age.
  for (int a = 0; a < kAges; ++a) {
    pop.numbers_at_age[0][a] = 1000.0;
  }

  AnnualTimeStep{}(data, life, fleet, 1000.0, pop);
  AnnualCatchAtStartOfYear{}(data, life, fleet, pop, fleet);

  constexpr double expected_y1_age1 = 1000.0;
  constexpr double expected_y1_age2 = 711.9177620339338;
  constexpr double expected_y1_age10 = 1226.9201954177395;
  constexpr double expected_catch_y0_age1 = 3.0479277436933629;
  constexpr double expected_catch_y1_age2 = 6.89337826858519;

  if (!nearly_equal(pop.numbers_at_age[1][0], expected_y1_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age1 got " << pop.numbers_at_age[1][0]
              << " expected " << expected_y1_age1 << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.numbers_at_age[1][1], expected_y1_age2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age2 got " << pop.numbers_at_age[1][1]
              << " expected " << expected_y1_age2
              << " diff " << (pop.numbers_at_age[1][1] - expected_y1_age2)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.numbers_at_age[1][9], expected_y1_age10)) {
    std::cerr << std::setprecision(17)
              << "FAIL: y1 age10 got " << pop.numbers_at_age[1][9]
              << " expected " << expected_y1_age10
              << " diff " << (pop.numbers_at_age[1][9] - expected_y1_age10)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.catch_numbers_at_age[0][0],
                    expected_catch_y0_age1)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch y0 age1 got "
              << fleet.catch_numbers_at_age[0][0]
              << " expected " << expected_catch_y0_age1
              << " diff "
              << (fleet.catch_numbers_at_age[0][0] - expected_catch_y0_age1)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fleet.catch_numbers_at_age[1][1],
                    expected_catch_y1_age2)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch y1 age2 got "
              << fleet.catch_numbers_at_age[1][1]
              << " expected " << expected_catch_y1_age2
              << " diff "
              << (fleet.catch_numbers_at_age[1][1] - expected_catch_y1_age2)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level12 annual timestep regression\n";
  return EXIT_SUCCESS;
}
