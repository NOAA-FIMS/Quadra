#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"

namespace bigeye_v2 {

//------------------------------------------------------------
// FishingMortality
//
// Purpose
// -------
// Computes fishing mortality-at-age and total mortality-at-age
// for one fleet.
//
// Consumes
// --------
// FleetParameters
// LifeHistoryState
// FleetState::selectivity_at_age
//
// Produces
// --------
// FleetState::f_at_age
// FleetState::z_at_age
//
// Notes
// -----
// Stateless.
// Owns no memory.
//------------------------------------------------------------
struct FishingMortality {
  template <typename T>
  void operator()(const FleetParameters<T> &p, const LifeHistoryState<T> &life,
                  FleetState<T> &fleet) const {
    for (int a = 0; a < kAges; ++a) {
      fleet.f_at_age[a] = p.fbar * fleet.selectivity_at_age[a];
      fleet.z_at_age[a] = life.m_at_age[a] + fleet.f_at_age[a];
    }
  }
};

} // namespace bigeye_v2
