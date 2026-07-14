#include "../common/joint_objective.hpp"

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

  data.observed_catch_biomass_by_year = {1326.1639007786976, 1350.0, 1300.0};

  data.observed_index_by_year = {116.46701723019194, 120.0, 110.0};

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
  p.catch_sigma = 0.1;
  p.q_index = 0.01;
  p.index_sigma = 0.2;

  BigeyeDerived<double> d;

  BigeyeJointObjective{}(data, p, d);

  constexpr double expected_catch_nll = 14.695989739827853;
  constexpr double expected_index_nll = 9.4692241297627753;
  constexpr double expected_agecomp_nll = 683.1930939334195;
  constexpr double expected_total_nll =
      expected_catch_nll + expected_index_nll + expected_agecomp_nll;

  if (!nearly_equal(d.catch_nll, expected_catch_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: catch_nll got " << d.catch_nll
              << " expected " << expected_catch_nll << " diff "
              << (d.catch_nll - expected_catch_nll) << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.index_nll, expected_index_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: index_nll got " << d.index_nll
              << " expected " << expected_index_nll << " diff "
              << (d.index_nll - expected_index_nll) << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.agecomp_nll, expected_agecomp_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: agecomp_nll got "
              << d.agecomp_nll << " expected " << expected_agecomp_nll
              << " diff " << (d.agecomp_nll - expected_agecomp_nll) << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(d.total_nll, expected_total_nll)) {
    std::cerr << std::setprecision(17) << "FAIL: total_nll got " << d.total_nll
              << " expected " << expected_total_nll << " diff "
              << (d.total_nll - expected_total_nll) << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 Level08 joint objective regression\n";
  return EXIT_SUCCESS;
}
