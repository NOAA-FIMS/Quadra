#include "../architecture/assessment/assessment_cycle.hpp"

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

  AssessmentParameters<double> parameters;
  parameters.life.log_m_young_offset = std::log(0.75);
  parameters.life.log_m_old_offset = std::log(0.65);

  parameters.populations.resize(1);
  parameters.populations[0].r0 = 1000.0;

  parameters.fleets.resize(1);
  parameters.fleets[0].sel_a50 = 5.0;
  parameters.fleets[0].sel_slope = 1.0;
  parameters.fleets[0].fbar = 0.2;
  parameters.fleets[0].q_index = 0.01;
  parameters.fleets[0].catch_sigma = 0.1;
  parameters.fleets[0].index_sigma = 0.2;

  AssessmentState<double> state;
  state.populations.resize(1);
  state.fleets.resize(1);

  state.populations[0].numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    state.populations[0].numbers_at_age[0][a] = 1000.0;
  }

  // First pass: generate synthetic observations from the CAA cycle.
  AssessmentCycle{}(data, parameters, state);

  // Use only years with positive biomass-index predictions.
  data.observed_catch_biomass_by_year = {
      state.fleets[0].total_catch_biomass_by_year[0],
      state.fleets[0].total_catch_biomass_by_year[1]};

  data.observed_index_by_year = {
      state.fleets[0].predicted_index_by_year[0],
      state.fleets[0].predicted_index_by_year[1]};

  data.catch_agecomp_sample_size = {100.0, 100.0};

  data.observed_catch_age_proportion = {
      state.fleets[0].predicted_catch_age_proportion[0],
      state.fleets[0].predicted_catch_age_proportion[1]};

  // Second pass: evaluate likelihood against self-generated observations.
  AssessmentState<double> fit_state;
  fit_state.populations.resize(1);
  fit_state.fleets.resize(1);

  fit_state.populations[0].numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    fit_state.populations[0].numbers_at_age[0][a] = 1000.0;
  }

  AssessmentCycle{}(data, parameters, fit_state);

  constexpr double expected_catch_nll = 13.868490024121582;
  constexpr double expected_index_nll = 10.080878113511577;
  constexpr double expected_agecomp_nll = 391.63707154283196;
  constexpr double expected_total_nll =
      expected_catch_nll + expected_index_nll + expected_agecomp_nll;

  if (!nearly_equal(fit_state.likelihood.catch_nll, expected_catch_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch_nll got " << fit_state.likelihood.catch_nll
              << " expected " << expected_catch_nll
              << " diff " << (fit_state.likelihood.catch_nll - expected_catch_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fit_state.likelihood.index_nll, expected_index_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: index_nll got " << fit_state.likelihood.index_nll
              << " expected " << expected_index_nll
              << " diff " << (fit_state.likelihood.index_nll - expected_index_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fit_state.likelihood.agecomp_nll, expected_agecomp_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: agecomp_nll got " << fit_state.likelihood.agecomp_nll
              << " expected " << expected_agecomp_nll
              << " diff "
              << (fit_state.likelihood.agecomp_nll - expected_agecomp_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(fit_state.likelihood.total_nll, expected_total_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: total_nll got " << fit_state.likelihood.total_nll
              << " expected " << expected_total_nll
              << " diff " << (fit_state.likelihood.total_nll - expected_total_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << std::setprecision(17);
  std::cout << "catch_nll," << fit_state.likelihood.catch_nll << "\n";
  std::cout << "index_nll," << fit_state.likelihood.index_nll << "\n";
  std::cout << "agecomp_nll," << fit_state.likelihood.agecomp_nll << "\n";
  std::cout << "total_nll," << fit_state.likelihood.total_nll << "\n";
  std::cout << "PASSED: Bigeye v2 CAA AssessmentCycle objective regression\n";

  return EXIT_SUCCESS;
}
