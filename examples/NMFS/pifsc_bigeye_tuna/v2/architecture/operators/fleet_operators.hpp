#pragma once
#include "../../common/bigeye_constants.hpp"
#include "../../common/derived.hpp"
#include "../../common/model_data.hpp"
#include "../parameters/fleet_parameters.hpp"
#include "../state/fleet_state.hpp"
#include "../state/population_state.hpp"
#include <cmath>

namespace bigeye_v2 {

struct FleetLogisticSelectivity {
  template <typename T>
  void operator()(const BigeyeModelData<T> &, const FleetParameters<T> &p,
                  FleetState<T> &fleet) const {
    for (int a = 0; a < kAges; ++a) {
      const T age = T(a + 1);
      fleet.selectivity_at_age[a] =
          T(1.0) / (T(1.0) + std::exp(-p.sel_slope * (age - p.sel_a50)));
    }
  }
};

struct FleetFishingMortality {
  template <typename T>
  void operator()(const BigeyeModelData<T> &, const FleetParameters<T> &p,
                  const BigeyeDerived<T> &life, FleetState<T> &fleet) const {
    for (int a = 0; a < kAges; ++a) {
      fleet.f_at_age[a] = p.fbar * fleet.selectivity_at_age[a];
      fleet.z_at_age[a] = life.m_at_age[a] + fleet.f_at_age[a];
    }
  }
};

} // namespace bigeye_v2
