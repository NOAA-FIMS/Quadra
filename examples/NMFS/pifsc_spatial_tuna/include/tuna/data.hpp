#ifndef QUADRA_MODEL_DATA_HPP
#define QUADRA_MODEL_DATA_HPP
#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace quadra {

struct TunaAssessmentData {
  int n_years_m = 0;
  int n_ages_m = 0;
  int n_fleets_m = 0;

  // Length n_ages_m
  std::vector<double> natural_mortality_at_age_m;
  std::vector<double> maturity_at_age_m;

  // Flattened as [year, age], length n_years_m * n_ages_m
  std::vector<double> weight_at_age_m;

  // Flattened as [fleet, year], length n_fleets_m * n_years_m
  std::vector<double> effort_m;
  std::vector<double> observed_index_m;

  // Flattened as [fleet, year, age], length n_fleets_m * n_years_m * n_ages_m
  std::vector<int> observed_catch_numbers_m;

  // Mid-year spawning fraction for spawning biomass calculation.
  double spawning_fraction_m = 0.5;

  size_t year_age_index(int year, int age) const {
    return static_cast<size_t>(year) * static_cast<size_t>(n_ages_m) +
           static_cast<size_t>(age);
  }

  size_t fleet_year_index(int fleet, int year) const {
    return static_cast<size_t>(fleet) * static_cast<size_t>(n_years_m) +
           static_cast<size_t>(year);
  }

  size_t fleet_year_age_index(int fleet, int year, int age) const {
    return (static_cast<size_t>(fleet) * static_cast<size_t>(n_years_m) +
            static_cast<size_t>(year)) *
               static_cast<size_t>(n_ages_m) +
           static_cast<size_t>(age);
  }

  void validate() const {
    if (n_years_m <= 1) {
      throw std::invalid_argument(
          "TunaAssessmentData: n_years_m must be at least 2");
    }

    if (n_ages_m <= 1) {
      throw std::invalid_argument(
          "TunaAssessmentData: n_ages_m must be at least 2");
    }

    if (n_fleets_m <= 0) {
      throw std::invalid_argument(
          "TunaAssessmentData: n_fleets_m must be positive");
    }

    const size_t expected_age = static_cast<size_t>(n_ages_m);
    const size_t expected_year_age =
        static_cast<size_t>(n_years_m) * static_cast<size_t>(n_ages_m);
    const size_t expected_fleet_year =
        static_cast<size_t>(n_fleets_m) * static_cast<size_t>(n_years_m);
    const size_t expected_fleet_year_age =
        expected_fleet_year * static_cast<size_t>(n_ages_m);

    if (natural_mortality_at_age_m.size() != expected_age) {
      throw std::invalid_argument(
          "TunaAssessmentData: natural_mortality_at_age_m has wrong length");
    }

    if (maturity_at_age_m.size() != expected_age) {
      throw std::invalid_argument(
          "TunaAssessmentData: maturity_at_age_m has wrong length");
    }

    if (weight_at_age_m.size() != expected_year_age) {
      throw std::invalid_argument(
          "TunaAssessmentData: weight_at_age_m has wrong length");
    }

    if (effort_m.size() != expected_fleet_year) {
      throw std::invalid_argument(
          "TunaAssessmentData: effort_m has wrong length");
    }

    if (observed_index_m.size() != expected_fleet_year) {
      throw std::invalid_argument(
          "TunaAssessmentData: observed_index_m has wrong length");
    }

    if (observed_catch_numbers_m.size() != expected_fleet_year_age) {
      throw std::invalid_argument(
          "TunaAssessmentData: observed_catch_numbers_m has wrong length");
    }

    if (spawning_fraction_m <= 0.0 || spawning_fraction_m >= 1.0) {
      throw std::invalid_argument(
          "TunaAssessmentData: spawning_fraction_m must be in (0,1)");
    }
  }
};

} // namespace quadra

#endif // QUADRA_MODEL_DATA_HPP
