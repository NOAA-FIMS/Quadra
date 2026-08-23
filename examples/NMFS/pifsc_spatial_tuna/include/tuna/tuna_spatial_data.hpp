#ifndef QUADRA_TUNA_SPATIAL_DATA_HPP
#define QUADRA_TUNA_SPATIAL_DATA_HPP
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace quadra {

struct TunaSpatialAssessmentData {
  int n_years_m = 0;
  int n_ages_m = 0;
  int n_fleets_m = 0;
  int n_regions_m = 0;
  int n_seasons_m = 4;

  // Length n_ages_m
  std::vector<double> natural_mortality_at_age_m;
  std::vector<double> maturity_at_age_m;

  // Flattened as [year, age], length n_years_m * n_ages_m
  std::vector<double> weight_at_age_m;

  // Flattened as [season, from_region, to_region],
  // length n_seasons_m * n_regions_m * n_regions_m
  std::vector<double> movement_matrix_m;

  // Length n_regions_m, sum should be 1.0
  std::vector<double> regional_recruit_proportions_m;

  // Flattened as [fleet, year, season, region],
  // length n_fleets_m * n_years_m * n_seasons_m * n_regions_m
  std::vector<double> effort_m;
  std::vector<double> observed_index_m;
  std::vector<double> observed_retained_biomass_m;
  std::vector<double> observed_discard_biomass_m;

  // Optional catch-conditioned removals, flattened as
  // [fleet, year, season, region]. catch_units_m has length n_fleets_m:
  // 1 = biomass, 2 = numbers. Empty vectors disable catch conditioning.
  std::vector<double> observed_total_catch_m;
  std::vector<int> catch_units_m;

  // Optional fixed availability surface flattened as [fleet, season, region,
  // age], length n_fleets_m * n_seasons_m * n_regions_m * n_ages_m.
  // If empty, the model assumes availability of 1.0.
  std::vector<double> availability_surface_m;

  // Flattened as [fleet, year, season, region, age],
  // length n_fleets_m * n_years_m * n_seasons_m * n_regions_m * n_ages_m
  std::vector<int> observed_catch_numbers_m;

  // Mid-year spawning fraction in annual spawning biomass calculation.
  double spawning_fraction_m = 0.5;

  size_t year_age_index(int year, int age) const {
    return static_cast<size_t>(year) * static_cast<size_t>(n_ages_m) +
           static_cast<size_t>(age);
  }

  size_t season_region_region_index(int season, int from_region,
                                    int to_region) const {
    return (static_cast<size_t>(season) * static_cast<size_t>(n_regions_m) +
            static_cast<size_t>(from_region)) *
               static_cast<size_t>(n_regions_m) +
           static_cast<size_t>(to_region);
  }

  size_t fleet_year_season_region_index(int fleet, int year, int season,
                                        int region) const {
    return ((static_cast<size_t>(fleet) * static_cast<size_t>(n_years_m) +
             static_cast<size_t>(year)) *
                static_cast<size_t>(n_seasons_m) +
            static_cast<size_t>(season)) *
               static_cast<size_t>(n_regions_m) +
           static_cast<size_t>(region);
  }

  size_t fleet_year_season_region_age_index(int fleet, int year, int season,
                                            int region, int age) const {
    return (fleet_year_season_region_index(fleet, year, season, region) *
            static_cast<size_t>(n_ages_m)) +
           static_cast<size_t>(age);
  }

  size_t fleet_season_region_age_index(int fleet, int season, int region,
                                       int age) const {
    return (((static_cast<size_t>(fleet) * static_cast<size_t>(n_seasons_m) +
              static_cast<size_t>(season)) *
                 static_cast<size_t>(n_regions_m) +
             static_cast<size_t>(region)) *
            static_cast<size_t>(n_ages_m)) +
           static_cast<size_t>(age);
  }

  size_t region_age_index(int region, int age) const {
    return static_cast<size_t>(region) * static_cast<size_t>(n_ages_m) +
           static_cast<size_t>(age);
  }

  void validate() const {
    const auto require_finite_range = [](const std::vector<double> &values,
                                         const std::string &name, double lower,
                                         double upper, bool lower_inclusive,
                                         bool upper_inclusive) {
      for (size_t i = 0; i < values.size(); ++i) {
        const double value = values[i];
        const bool lower_ok = lower_inclusive ? value >= lower : value > lower;
        const bool upper_ok = upper_inclusive ? value <= upper : value < upper;
        if (!std::isfinite(value) || !lower_ok || !upper_ok) {
          throw std::invalid_argument("TunaSpatialAssessmentData: " + name +
                                      "[" + std::to_string(i) +
                                      "] is outside its valid range");
        }
      }
    };
    const double infinity = std::numeric_limits<double>::infinity();
    if (n_years_m <= 1) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: n_years_m must be at least 2");
    }

    if (n_ages_m <= 1) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: n_ages_m must be at least 2");
    }

    if (n_fleets_m <= 0) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: n_fleets_m must be positive");
    }

    if (n_regions_m <= 0) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: n_regions_m must be positive");
    }

    if (n_seasons_m <= 0) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: n_seasons_m must be positive");
    }

    if (spawning_fraction_m <= 0.0 || spawning_fraction_m >= 1.0) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: spawning_fraction_m must be in (0,1)");
    }

    const size_t expected_age = static_cast<size_t>(n_ages_m);
    const size_t expected_year_age =
        static_cast<size_t>(n_years_m) * static_cast<size_t>(n_ages_m);
    const size_t expected_movement = static_cast<size_t>(n_seasons_m) *
                                     static_cast<size_t>(n_regions_m) *
                                     static_cast<size_t>(n_regions_m);
    const size_t expected_fysr =
        static_cast<size_t>(n_fleets_m) * static_cast<size_t>(n_years_m) *
        static_cast<size_t>(n_seasons_m) * static_cast<size_t>(n_regions_m);
    const size_t expected_fysra = expected_fysr * static_cast<size_t>(n_ages_m);
    const size_t expected_fsra =
        static_cast<size_t>(n_fleets_m) * static_cast<size_t>(n_seasons_m) *
        static_cast<size_t>(n_regions_m) * static_cast<size_t>(n_ages_m);

    if (natural_mortality_at_age_m.size() != expected_age) {
      throw std::invalid_argument("TunaSpatialAssessmentData: "
                                  "natural_mortality_at_age_m has wrong "
                                  "length");
    }

    if (maturity_at_age_m.size() != expected_age) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: maturity_at_age_m has wrong length");
    }

    if (weight_at_age_m.size() != expected_year_age) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: weight_at_age_m has wrong length");
    }

    if (movement_matrix_m.size() != expected_movement) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: movement_matrix_m has wrong length");
    }

    if (regional_recruit_proportions_m.size() !=
        static_cast<size_t>(n_regions_m)) {
      throw std::invalid_argument("TunaSpatialAssessmentData: "
                                  "regional_recruit_proportions_m has wrong "
                                  "length");
    }

    if (effort_m.size() != expected_fysr) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: effort_m has wrong length");
    }

    if (observed_index_m.size() != expected_fysr) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: observed_index_m has wrong length");
    }

    if (observed_retained_biomass_m.size() != expected_fysr) {
      throw std::invalid_argument("TunaSpatialAssessmentData: "
                                  "observed_retained_biomass_m has wrong "
                                  "length");
    }

    if (observed_discard_biomass_m.size() != expected_fysr) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: observed_discard_biomass_m has wrong "
          "length");
    }

    if (!observed_total_catch_m.empty() &&
        observed_total_catch_m.size() != expected_fysr) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: observed_total_catch_m has wrong length");
    }
    if (!catch_units_m.empty() &&
        catch_units_m.size() != static_cast<size_t>(n_fleets_m)) {
      throw std::invalid_argument(
          "TunaSpatialAssessmentData: catch_units_m has wrong length");
    }
    if (observed_total_catch_m.empty() != catch_units_m.empty()) {
      throw std::invalid_argument("TunaSpatialAssessmentData: conditioned "
                                  "catch and units must both be supplied");
    }
    for (int unit : catch_units_m) {
      if (unit != 1 && unit != 2)
        throw std::invalid_argument("TunaSpatialAssessmentData: catch units "
                                    "must be 1 (biomass) or 2 (numbers)");
    }

    if (observed_catch_numbers_m.size() != expected_fysra) {
      throw std::invalid_argument("TunaSpatialAssessmentData: "
                                  "observed_catch_numbers_m has wrong length");
    }

    if (!availability_surface_m.empty() &&
        availability_surface_m.size() != expected_fsra) {
      throw std::invalid_argument("TunaSpatialAssessmentData: "
                                  "availability_surface_m has wrong length");
    }

    require_finite_range(natural_mortality_at_age_m,
                         "natural_mortality_at_age_m", 0.0, infinity, true,
                         true);
    require_finite_range(maturity_at_age_m, "maturity_at_age_m", 0.0, 1.0, true,
                         true);
    require_finite_range(weight_at_age_m, "weight_at_age_m", 0.0, infinity,
                         false, true);
    require_finite_range(effort_m, "effort_m", 0.0, infinity, true, true);
    require_finite_range(observed_index_m, "observed_index_m", 0.0, infinity,
                         true, true);
    require_finite_range(observed_retained_biomass_m,
                         "observed_retained_biomass_m", 0.0, infinity, true,
                         true);
    require_finite_range(observed_discard_biomass_m,
                         "observed_discard_biomass_m", 0.0, infinity, true,
                         true);
    if (!observed_total_catch_m.empty()) {
      require_finite_range(observed_total_catch_m, "observed_total_catch_m",
                           0.0, infinity, true, true);
    }
    if (!availability_surface_m.empty()) {
      require_finite_range(availability_surface_m, "availability_surface_m",
                           0.0, infinity, true, true);
    }
    for (size_t i = 0; i < observed_catch_numbers_m.size(); ++i) {
      if (observed_catch_numbers_m[i] < 0) {
        throw std::invalid_argument(
            "TunaSpatialAssessmentData: observed_catch_numbers_m[" +
            std::to_string(i) + "] must be non-negative");
      }
    }

    double recruit_sum = 0.0;
    for (double p : regional_recruit_proportions_m) {
      if (!std::isfinite(p) || p < 0.0) {
        throw std::invalid_argument("TunaSpatialAssessmentData: "
                                    "regional_recruit_proportions_m must be "
                                    "non-negative");
      }
      recruit_sum += p;
    }

    if (std::fabs(recruit_sum - 1.0) > 1e-8) {
      throw std::invalid_argument("TunaSpatialAssessmentData: "
                                  "regional_recruit_proportions_m must sum "
                                  "to 1.0");
    }

    for (int s = 0; s < n_seasons_m; ++s) {
      for (int from = 0; from < n_regions_m; ++from) {
        double row_sum = 0.0;
        for (int to = 0; to < n_regions_m; ++to) {
          const double p =
              movement_matrix_m[season_region_region_index(s, from, to)];
          if (!std::isfinite(p) || p < 0.0) {
            throw std::invalid_argument(
                "TunaSpatialAssessmentData: movement probabilities must be "
                "non-negative");
          }
          row_sum += p;
        }

        if (std::fabs(row_sum - 1.0) > 1e-8) {
          throw std::invalid_argument(
              "TunaSpatialAssessmentData: movement rows must sum to 1.0");
        }
      }
    }
  }

  TunaSpatialAssessmentData peel_years(int n_years_keep) const {
    if (n_years_keep < 2 || n_years_keep > n_years_m) {
      throw std::invalid_argument("TunaSpatialAssessmentData::peel_years: "
                                  "n_years_keep must be in [2, n_years_m]");
    }

    TunaSpatialAssessmentData out = *this;
    out.n_years_m = n_years_keep;

    out.weight_at_age_m.assign(static_cast<size_t>(out.n_years_m) *
                                   static_cast<size_t>(out.n_ages_m),
                               0.0);

    for (int y = 0; y < out.n_years_m; ++y) {
      for (int a = 0; a < out.n_ages_m; ++a) {
        out.weight_at_age_m[out.year_age_index(y, a)] =
            weight_at_age_m[year_age_index(y, a)];
      }
    }

    const size_t n_fysr = static_cast<size_t>(out.n_fleets_m) *
                          static_cast<size_t>(out.n_years_m) *
                          static_cast<size_t>(out.n_seasons_m) *
                          static_cast<size_t>(out.n_regions_m);
    out.effort_m.assign(n_fysr, 0.0);
    out.observed_index_m.assign(n_fysr, 0.0);
    out.observed_retained_biomass_m.assign(n_fysr, 0.0);
    out.observed_discard_biomass_m.assign(n_fysr, 0.0);
    if (!observed_total_catch_m.empty())
      out.observed_total_catch_m.assign(n_fysr, 0.0);

    for (int f = 0; f < out.n_fleets_m; ++f) {
      for (int y = 0; y < out.n_years_m; ++y) {
        for (int s = 0; s < out.n_seasons_m; ++s) {
          for (int r = 0; r < out.n_regions_m; ++r) {
            const size_t src = fleet_year_season_region_index(f, y, s, r);
            const size_t dst = out.fleet_year_season_region_index(f, y, s, r);
            out.effort_m[dst] = effort_m[src];
            out.observed_index_m[dst] = observed_index_m[src];
            out.observed_retained_biomass_m[dst] =
                observed_retained_biomass_m[src];
            out.observed_discard_biomass_m[dst] =
                observed_discard_biomass_m[src];
            if (!observed_total_catch_m.empty())
              out.observed_total_catch_m[dst] = observed_total_catch_m[src];
          }
        }
      }
    }

    const size_t n_fysra = n_fysr * static_cast<size_t>(out.n_ages_m);
    out.observed_catch_numbers_m.assign(n_fysra, 0);

    for (int f = 0; f < out.n_fleets_m; ++f) {
      for (int y = 0; y < out.n_years_m; ++y) {
        for (int s = 0; s < out.n_seasons_m; ++s) {
          for (int r = 0; r < out.n_regions_m; ++r) {
            for (int a = 0; a < out.n_ages_m; ++a) {
              const size_t src =
                  fleet_year_season_region_age_index(f, y, s, r, a);
              const size_t dst =
                  out.fleet_year_season_region_age_index(f, y, s, r, a);
              out.observed_catch_numbers_m[dst] = observed_catch_numbers_m[src];
            }
          }
        }
      }
    }

    out.validate();
    return out;
  }
};

} // namespace quadra

#endif // QUADRA_TUNA_SPATIAL_DATA_HPP
