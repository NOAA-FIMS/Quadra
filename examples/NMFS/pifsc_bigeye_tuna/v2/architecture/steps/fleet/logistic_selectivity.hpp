#pragma once

#include "../../../common/bigeye_constants.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../../common/model_data.hpp"

#include <cmath>

namespace bigeye_v2 {

// Logistic selectivity.
// Computes vulnerability-at-age for one fleet.
struct LogisticSelectivity {
  template <typename T>
  void operator()(const BigeyeModelData<T> &,
                  const FleetParameters<T> &p,
                  FleetState<T> &fleet) const {
    for (int a = 0; a < kAges; ++a) {
      const T age = T(a + 1);
      fleet.selectivity_at_age[a] =
          T(1.0) / (T(1.0) + std::exp(-p.sel_slope * (age - p.sel_a50)));
    }
  }
};

} // namespace bigeye_v2
