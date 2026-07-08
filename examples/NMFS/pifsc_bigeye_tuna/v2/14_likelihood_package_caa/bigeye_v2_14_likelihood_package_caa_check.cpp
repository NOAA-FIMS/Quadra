#include "../architecture/packages/likelihood/likelihood_package.hpp"

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

  data.observed_catch_biomass_by_year = {100.0, 120.0, 90.0};
  data.observed_index_by_year = {10.0, 12.0, 9.0};
  data.catch_agecomp_sample_size = {100.0, 100.0, 100.0};

  data.observed_catch_age_proportion = {
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06}};

  FleetParameters<double> parameters;
  parameters.catch_sigma = 0.1;
  parameters.index_sigma = 0.2;

  FleetState<double> fleet;
  fleet.total_catch_biomass_by_year = {100.0, 120.0, 90.0};
  fleet.predicted_index_by_year = {10.0, 12.0, 9.0};

  fleet.predicted_catch_age_proportion.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int y = 0; y < data.n_years; ++y) {
    fleet.predicted_catch_age_proportion[y] =
        std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                  0.16, 0.17, 0.14, 0.10, 0.06};
  }

  LikelihoodState<double> likelihood;

  LikelihoodPackage{}(data, parameters, fleet, likelihood);

  constexpr double expected_catch_nll = 6.9847163201182667;
  constexpr double expected_index_nll = 2.1564025828159652;
  constexpr double expected_agecomp_nll = 650.2267963186614;
  constexpr double expected_total_nll =
      expected_catch_nll + expected_index_nll + expected_agecomp_nll;

  if (!nearly_equal(likelihood.catch_nll, expected_catch_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch_nll got " << likelihood.catch_nll
              << " expected " << expected_catch_nll
              << " diff " << (likelihood.catch_nll - expected_catch_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(likelihood.index_nll, expected_index_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: index_nll got " << likelihood.index_nll
              << " expected " << expected_index_nll
              << " diff " << (likelihood.index_nll - expected_index_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(likelihood.agecomp_nll, expected_agecomp_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: agecomp_nll got " << likelihood.agecomp_nll
              << " expected " << expected_agecomp_nll
              << " diff " << (likelihood.agecomp_nll - expected_agecomp_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(likelihood.total_nll, expected_total_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: total_nll got " << likelihood.total_nll
              << " expected " << expected_total_nll
              << " diff " << (likelihood.total_nll - expected_total_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA LikelihoodPackage regression\n";
  return EXIT_SUCCESS;
}
