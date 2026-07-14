#pragma once

#include "../../../common/bigeye_constants.hpp"
#include "../../../common/model_data.hpp"
#include "../../state/fleet_state.hpp"

#include <array>

namespace bigeye_v2 {

// Predicts catch age composition from fleet catch numbers-at-age.
struct CatchAgeCompositionPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data, FleetState<T> &fleet) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    fleet.predicted_catch_age_proportion.assign(ny, std::array<T, kAges>{});

    for (std::size_t y = 0; y < ny; ++y) {
      T total = T(0.0);

      for (int a = 0; a < kAges; ++a) {
        total += fleet.catch_numbers_at_age[y][a];
      }

      for (int a = 0; a < kAges; ++a) {
        fleet.predicted_catch_age_proportion[y][a] =
            fleet.catch_numbers_at_age[y][a] / total;
      }
    }
  }
};

} // namespace bigeye_v2
