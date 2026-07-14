#include "../architecture/packages/population/population_package.hpp"
#include "../architecture/parameters/fleet_parameters.hpp"
#include "../architecture/parameters/life_history_parameters.hpp"
#include "../architecture/state/fleet_state.hpp"
#include "../architecture/state/life_history_state.hpp"
#include "../architecture/steps/fleet/fishing_mortality.hpp"
#include "../architecture/steps/fleet/logistic_selectivity.hpp"
#include "../architecture/steps/life_history/life_history.hpp"
#include "../common/model_data.hpp"

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

  LifeHistoryParameters<double> lp;
  lp.log_m_young_offset = std::log(0.75);
  lp.log_m_old_offset = std::log(0.65);

  LifeHistoryState<double> life;
  BigeyeLifeHistory{}(data, lp, life);

  FleetParameters<double> fp;
  fp.sel_a50 = 5.0;
  fp.sel_slope = 1.0;
  fp.fbar = 0.2;

  FleetState<double> fleet;
  LogisticSelectivity{}(data, fp, fleet);
  FishingMortality{}(fp, life, fleet);

  PopulationParameters<double> pp;
  pp.r0 = 1000.0;

  PopulationState<double> pop;

  // Initial numbers for year 0. Package keeps age-1 recruitment and
  // advances survivors into later years.
  pop.numbers_at_age.assign(static_cast<std::size_t>(data.n_years),
                            std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    pop.numbers_at_age[0][a] = 1000.0;
  }

  PopulationPackage{}(data, pp, life, fleet, pop);

  constexpr double expected_y0_age1 = 1000.0;
  constexpr double expected_y1_age1 = 1000.0;
  constexpr double expected_y1_age2 = 710.98976678322833;
  constexpr double expected_y1_age10 = 1225.2142471639108;
  constexpr double expected_ssb_y0 = 91000.0;

  if (!nearly_equal(pop.numbers_at_age[0][0], expected_y0_age1)) {
    std::cerr << "FAIL: y0 age1\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.numbers_at_age[1][0], expected_y1_age1)) {
    std::cerr << "FAIL: y1 age1\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.numbers_at_age[1][1], expected_y1_age2)) {
    std::cerr << std::setprecision(17) << "FAIL: y1 age2 got "
              << pop.numbers_at_age[1][1] << " expected " << expected_y1_age2
              << " diff " << (pop.numbers_at_age[1][1] - expected_y1_age2)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.numbers_at_age[1][9], expected_y1_age10)) {
    std::cerr << std::setprecision(17) << "FAIL: y1 age10 got "
              << pop.numbers_at_age[1][9] << " expected " << expected_y1_age10
              << " diff " << (pop.numbers_at_age[1][9] - expected_y1_age10)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(pop.spawning_biomass_by_year[0], expected_ssb_y0)) {
    std::cerr << std::setprecision(17) << "FAIL: ssb y0 got "
              << pop.spawning_biomass_by_year[0] << " expected "
              << expected_ssb_y0 << " diff "
              << (pop.spawning_biomass_by_year[0] - expected_ssb_y0) << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA PopulationPackage regression\n";
  return EXIT_SUCCESS;
}
