#pragma once

#include "../../common/bigeye_constants.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T> struct FleetState {
  std::array<T, kAges> selectivity_at_age{};
  std::array<T, kAges> f_at_age{};
  std::array<T, kAges> z_at_age{};

  std::vector<std::array<T, kAges>> catch_numbers_at_age{};
  std::vector<std::array<T, kAges>> catch_biomass_at_age{};
  std::vector<T> total_catch_biomass_by_year{};
  std::vector<T> predicted_index_by_year{};
  std::vector<std::array<T, kAges>> predicted_catch_age_proportion{};
};

} // namespace bigeye_v2
