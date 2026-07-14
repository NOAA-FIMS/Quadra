#include "../common/catch.hpp"
#include "../common/index.hpp"
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
  data.observed_index_by_year = {116.46701723019194, 120.0, 110.0};

  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);
  p.r0 = 1000.0;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.fbar = 0.2;
  p.q_index = 0.01;
  p.index_sigma = 0.2;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, p, d);
  FixedRecruitment{}(data, p, d);
  LogisticSelectivity{}(data, p, d);
  FishingMortality{}(data, p, d);
  UnfishedPopulation{}(data, p, d);
  BiomassIndexPrediction{}(data, p, d);
  LognormalIndexLikelihood{}(data, p, d);

  constexpr double expected_pred_index = 116.46701723019194;
  constexpr double expected_index_nll = 9.4692241297627753;

  if (!nearly_equal(d.predicted_index_by_year[0], expected_pred_index)) {
    std::cerr << std::setprecision(17) << "FAIL: predicted_index got "
              << d.predicted_index_by_year[0] << " expected "
              << expected_pred_index << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.index_nll, expected_index_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: index_nll got " << d.index_nll
              << " expected " << expected_index_nll << " diff "
              << (d.index_nll - expected_index_nll) << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level06 index likelihood regression\n";
  return EXIT_SUCCESS;
}
