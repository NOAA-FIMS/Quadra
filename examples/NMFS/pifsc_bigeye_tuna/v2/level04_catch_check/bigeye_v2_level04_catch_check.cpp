#include "../common/catch.hpp"
#include "../common/life_history.hpp"
#include "../common/mortality.hpp"
#include "../common/population.hpp"
#include "../common/recruitment.hpp"
#include "../common/selectivity.hpp"

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
  p.r0 = 1000.0;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.fbar = 0.2;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  LogisticSelectivity{}(data, p, d);
  FishingMortality{}(data, p, d);
  UnfishedPopulation{}(data, p, d);
  BaranovCatch{}(data, p, d);

  constexpr double expected_catch_n_age_1 = 3.0479277436933629;
  constexpr double expected_catch_n_age_5 = 17.818597787065155;
  constexpr double expected_catch_n_age_10 = 20.708228383041849;
  constexpr double expected_total_catch_biomass = 1326.1639007786976;

  for (int y = 0; y < data.n_years; ++y) {
    if (!nearly_equal(d.catch_numbers_at_age[y][0], expected_catch_n_age_1)) {
      std::cerr << std::setprecision(17)
                << "FAIL: catch_numbers_at_age[" << y << "][0] got "
                << d.catch_numbers_at_age[y][0]
                << " expected " << expected_catch_n_age_1
                << " diff "
                << (d.catch_numbers_at_age[y][0] - expected_catch_n_age_1)
                << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.catch_numbers_at_age[y][4], expected_catch_n_age_5)) {
      std::cerr << std::setprecision(17)
                << "FAIL: catch_numbers_at_age[" << y << "][4] got "
                << d.catch_numbers_at_age[y][4]
                << " expected " << expected_catch_n_age_5
                << " diff "
                << (d.catch_numbers_at_age[y][4] - expected_catch_n_age_5)
                << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.catch_numbers_at_age[y][9], expected_catch_n_age_10)) {
      std::cerr << std::setprecision(17)
                << "FAIL: catch_numbers_at_age[" << y << "][9] got "
                << d.catch_numbers_at_age[y][9]
                << " expected " << expected_catch_n_age_10
                << " diff "
                << (d.catch_numbers_at_age[y][9] - expected_catch_n_age_10)
                << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.total_catch_biomass_by_year[y],
                      expected_total_catch_biomass)) {
      std::cerr << std::setprecision(17)
                << "FAIL: total_catch_biomass_by_year[" << y << "] got "
                << d.total_catch_biomass_by_year[y]
                << " expected " << expected_total_catch_biomass
                << " diff "
                << (d.total_catch_biomass_by_year[y] -
                    expected_total_catch_biomass)
                << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 Level04 catch regression\n";
  return EXIT_SUCCESS;
}
