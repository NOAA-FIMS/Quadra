#pragma once

#include "bigeye_constants.hpp"

#include <array>
#include <vector>

namespace bigeye_v2 {

template <typename T> struct BigeyeModelData {
  int n_years = 1;
  std::vector<T> observed_catch_biomass_by_year{};
  std::vector<T> observed_index_by_year{};
  std::vector<std::array<T, kAges>> observed_catch_age_proportion{};
  std::vector<T> catch_agecomp_sample_size{};
};

} // namespace bigeye_v2
