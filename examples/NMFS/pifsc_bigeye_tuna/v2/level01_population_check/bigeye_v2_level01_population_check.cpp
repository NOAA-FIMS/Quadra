#include "../common/life_history.hpp"
#include "../common/population.hpp"
#include "../common/recruitment.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {

bool nearly_equal(double a, double b, double tol = 1.0e-5) {
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

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  UnfishedPopulation{}(data, p, d);

  constexpr double expected_recruit = 1000.0;
  constexpr double expected_n_age_1 = 1000.0;
  constexpr double expected_n_age_2 = 713.55197470650239;
  constexpr double expected_n_age_10 = 131.92500084784737;
  constexpr double expected_ssb = 11646.701723019194;

  for (int y = 0; y < data.n_years; ++y) {
    if (!nearly_equal(d.recruits_by_year[y], expected_recruit)) {
      std::cerr << "FAIL: recruits_by_year[" << y << "]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.numbers_at_age[y][0], expected_n_age_1)) {
      std::cerr << "FAIL: numbers_at_age[" << y << "][0]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.numbers_at_age[y][1], expected_n_age_2)) {
      std::cerr << std::setprecision(17) << "FAIL: numbers_at_age[" << y
                << "][1] got " << d.numbers_at_age[y][1] << " expected "
                << expected_n_age_2 << " diff "
                << (d.numbers_at_age[y][1] - expected_n_age_2) << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.numbers_at_age[y][9], expected_n_age_10)) {
      std::cerr << std::setprecision(17) << "FAIL: numbers_at_age[" << y
                << "][9] got " << d.numbers_at_age[y][9] << " expected "
                << expected_n_age_10 << " diff "
                << (d.numbers_at_age[y][9] - expected_n_age_10) << "\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.spawning_biomass_by_year[y], expected_ssb)) {
      std::cerr << std::setprecision(17) << "FAIL: spawning_biomass_by_year["
                << y << "] got " << d.spawning_biomass_by_year[y]
                << " expected " << expected_ssb << " diff "
                << (d.spawning_biomass_by_year[y] - expected_ssb) << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 Level01 population regression\n";
  return EXIT_SUCCESS;
}
