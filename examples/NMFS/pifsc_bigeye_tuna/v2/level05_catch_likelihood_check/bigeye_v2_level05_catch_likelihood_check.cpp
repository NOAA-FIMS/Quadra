#include "../common/catch.hpp"
#include "../common/life_history.hpp"
#include "../common/likelihood.hpp"
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
  data.observed_catch_biomass_by_year = {1326.1639007786976, 1350.0, 1300.0};

  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);
  p.r0 = 1000.0;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.fbar = 0.2;
  p.catch_sigma = 0.1;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  LogisticSelectivity{}(data, p, d);
  FishingMortality{}(data, p, d);
  UnfishedPopulation{}(data, p, d);
  BaranovCatch{}(data, p, d);
  LognormalCatchLikelihood{}(data, p, d);

  constexpr double expected_catch_nll = 14.695989739827853;
  constexpr double expected_total_nll = 14.695989739827853;

  if (!nearly_equal(d.catch_nll, expected_catch_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: catch_nll got " << d.catch_nll
              << " expected " << expected_catch_nll << " diff "
              << (d.catch_nll - expected_catch_nll) << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.total_nll, expected_total_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: total_nll got " << d.total_nll
              << " expected " << expected_total_nll << " diff "
              << (d.total_nll - expected_total_nll) << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level05 catch likelihood regression\n";
  return EXIT_SUCCESS;
}
