#pragma once

#include "../../common/bigeye_constants.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T> struct PopulationState {
  std::vector<T> recruits_by_year{};
  std::vector<std::array<T, kAges>> numbers_at_age{};
  std::vector<std::array<T, kAges>> survivors_at_age{};
  std::vector<T> spawning_biomass_by_year{};
};

} // namespace bigeye_v2
