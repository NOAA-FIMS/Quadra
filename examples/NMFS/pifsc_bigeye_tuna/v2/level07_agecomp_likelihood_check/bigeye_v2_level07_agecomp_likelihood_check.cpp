#include "../common/agecomp.hpp"
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
  data.catch_agecomp_sample_size = {100.0, 100.0, 100.0};

  data.observed_catch_age_proportion = {
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14, 0.16, 0.17, 0.14,
                                0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14, 0.16, 0.17, 0.14,
                                0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14, 0.16, 0.17, 0.14,
                                0.10, 0.06}};

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
  CatchAgeCompositionPrediction{}(data, p, d);
  MultinomialAgeCompLikelihood{}(data, p, d);

  constexpr double expected_prop_age_1 = 0.025946689511201496;
  constexpr double expected_prop_age_10 = 0.17628697835557786;
  constexpr double expected_agecomp_nll = 683.1930939334195;

  if (!nearly_equal(d.predicted_catch_age_proportion[0][0],
                    expected_prop_age_1)) {
    std::cerr << std::setprecision(17) << "FAIL: prop age 1 got "
              << d.predicted_catch_age_proportion[0][0] << " expected "
              << expected_prop_age_1 << " diff "
              << (d.predicted_catch_age_proportion[0][0] - expected_prop_age_1)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.predicted_catch_age_proportion[0][9],
                    expected_prop_age_10)) {
    std::cerr << std::setprecision(17) << "FAIL: prop age 10 got "
              << d.predicted_catch_age_proportion[0][9] << " expected "
              << expected_prop_age_10 << " diff "
              << (d.predicted_catch_age_proportion[0][9] - expected_prop_age_10)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.agecomp_nll, expected_agecomp_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: agecomp_nll got "
              << d.agecomp_nll << " expected " << expected_agecomp_nll
              << " diff " << (d.agecomp_nll - expected_agecomp_nll) << "\n";
    return EXIT_FAILURE;
  }

  std::cout
      << "PASSED: Bigeye v2 Level07 age composition likelihood regression\n";
  return EXIT_SUCCESS;
}
