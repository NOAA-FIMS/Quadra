#include "../common/agecomp.hpp"
#include "../common/catch.hpp"
#include "../common/index.hpp"
#include "../common/joint_objective.hpp"
#include "../common/life_history.hpp"
#include "../common/mortality.hpp"
#include "../common/population.hpp"
#include "../common/recruitment.hpp"
#include "../common/selectivity.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace bigeye_v2;

namespace {

BigeyeModelParameters<double> base_params() {
  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);
  p.r0 = 1000.0;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;
  p.catch_sigma = 0.1;
  p.index_sigma = 0.2;
  return p;
}

BigeyeModelData<double> make_synthetic_data() {
  BigeyeModelData<double> data;
  data.n_years = 3;
  data.catch_agecomp_sample_size = {100.0, 100.0, 100.0};

  auto true_p = base_params();
  true_p.fbar = 0.2;
  true_p.q_index = 0.01;

  BigeyeDerived<double> d;

  BigeyeLifeHistory{}(data, true_p, d);
  FixedRecruitment{}(data, true_p, d);
  LogisticSelectivity{}(data, true_p, d);
  FishingMortality{}(data, true_p, d);
  UnfishedPopulation{}(data, true_p, d);
  BaranovCatch{}(data, true_p, d);
  BiomassIndexPrediction{}(data, true_p, d);
  CatchAgeCompositionPrediction{}(data, true_p, d);

  data.observed_catch_biomass_by_year = d.total_catch_biomass_by_year;
  data.observed_index_by_year = d.predicted_index_by_year;
  data.observed_catch_age_proportion = d.predicted_catch_age_proportion;

  return data;
}

double evaluate(const BigeyeModelData<double> &data, double log_fbar,
                double log_q) {
  auto p = base_params();
  p.fbar = std::exp(log_fbar);
  p.q_index = std::exp(log_q);

  BigeyeDerived<double> d;
  BigeyeJointObjective{}(data, p, d);
  return d.total_nll;
}

bool nearly_equal(double a, double b, double tol) {
  return std::abs(a - b) <= tol;
}

} // namespace

int main() {
  const auto data = make_synthetic_data();

  double log_fbar = std::log(0.08);
  double log_q = std::log(0.03);

  double best = evaluate(data, log_fbar, log_q);
  double step = 0.5;

  for (int iter = 0; iter < 200 && step > 1.0e-10; ++iter) {
    bool improved = false;

    const double candidates[4][2] = {
        {log_fbar + step, log_q},
        {log_fbar - step, log_q},
        {log_fbar, log_q + step},
        {log_fbar, log_q - step},
    };

    for (const auto &c : candidates) {
      const double fx = evaluate(data, c[0], c[1]);
      if (fx < best) {
        best = fx;
        log_fbar = c[0];
        log_q = c[1];
        improved = true;
      }
    }

    if (!improved) {
      step *= 0.5;
    }
  }

  const double fbar_hat = std::exp(log_fbar);
  const double q_hat = std::exp(log_q);

  std::cout << std::setprecision(17);
  std::cout << "fbar_hat," << fbar_hat << "\n";
  std::cout << "q_hat," << q_hat << "\n";
  std::cout << "objective," << best << "\n";

  if (!nearly_equal(fbar_hat, 0.2, 1.0e-3)) {
    std::cerr << "FAIL: fbar recovery\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(q_hat, 0.01, 1.0e-4)) {
    std::cerr << "FAIL: q recovery\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level09 optimizer recovery regression\n";
  return EXIT_SUCCESS;
}
