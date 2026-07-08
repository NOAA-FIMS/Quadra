#pragma once

#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"
#include "../../../common/model_data.hpp"

#include <array>
#include <cmath>

namespace bigeye_v2 {

//------------------------------------------------------------
// BaranovCatch
//
// Purpose
// -------
// Computes catch-at-age and catch biomass for one fleet.
//
// Consumes
// --------
// BigeyeModelData
// LifeHistoryState
// PopulationState::numbers_at_age
// FleetState::f_at_age
// FleetState::z_at_age
//
// Produces
// --------
// FleetState::catch_numbers_at_age
// FleetState::catch_biomass_at_age
// FleetState::total_catch_biomass_by_year
//
// Notes
// -----
// Stateless.
// Owns no memory beyond resizing fleet-owned state.
//------------------------------------------------------------
struct BaranovCatch {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const LifeHistoryState<T> &life,
                  const PopulationState<T> &population,
                  FleetState<T> &fleet) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    fleet.catch_numbers_at_age.assign(ny, std::array<T, kAges>{});
    fleet.catch_biomass_at_age.assign(ny, std::array<T, kAges>{});
    fleet.total_catch_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        const T z = fleet.z_at_age[a];

        const T cn =
            population.numbers_at_age[y][a] *
            (fleet.f_at_age[a] / z) *
            (T(1.0) - std::exp(-z));

        fleet.catch_numbers_at_age[y][a] = cn;
        fleet.catch_biomass_at_age[y][a] = cn * life.weight_at_age[a];
        fleet.total_catch_biomass_by_year[y] += fleet.catch_biomass_at_age[y][a];
      }
    }
  }
};

} // namespace bigeye_v2
