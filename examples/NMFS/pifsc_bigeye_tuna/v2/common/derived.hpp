#pragma once

#include "bigeye_constants.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T> struct BigeyeDerived {
  std::array<T, kAges> m_at_age{};
  std::array<T, kAges> weight_at_age{};
  std::array<T, kAges> maturity_at_age{};
  std::array<T, kAges> selectivity_at_age{};
  std::array<T, kAges> f_at_age{};
  std::array<T, kAges> z_at_age{};

  std::vector<T> recruits_by_year{};
  std::vector<std::array<T, kAges>> numbers_at_age{};
  std::vector<T> spawning_biomass_by_year{};

  std::vector<std::array<T, kAges>> catch_numbers_at_age{};
  std::vector<std::array<T, kAges>> catch_biomass_at_age{};
  std::vector<T> total_catch_biomass_by_year{};
  std::vector<T> predicted_index_by_year{};
  std::vector<std::array<T, kAges>> predicted_catch_age_proportion{};

  T catch_nll = T(0.0);
  T index_nll = T(0.0);
  T agecomp_nll = T(0.0);
  T total_nll = T(0.0);
};

} // namespace bigeye_v2
