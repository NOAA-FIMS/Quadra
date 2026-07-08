#pragma once

#include "../../../common/bigeye_constants.hpp"
#include "../../parameters/life_history_parameters.hpp"
#include "../../state/life_history_state.hpp"
#include "../../../common/model_data.hpp"

#include <cmath>

namespace bigeye_v2 {

// Bigeye life history.
// Computes natural mortality, weight, and maturity at age.
struct BigeyeLifeHistory {
  template <typename T>
  void operator()(const BigeyeModelData<T> &,
                  const LifeHistoryParameters<T> &p,
                  LifeHistoryState<T> &life) const {
    const T adult_m = T(0.45);
    const T m_young = adult_m * std::exp(p.log_m_young_offset);
    const T m_old = adult_m * std::exp(p.log_m_old_offset);

    for (int a = 0; a < kAges; ++a) {
      const int age = a + 1;

      if (age <= 3) {
        life.m_at_age[a] = m_young;
      } else if (age >= 8) {
        life.m_at_age[a] = m_old;
      } else {
        life.m_at_age[a] = adult_m;
      }

      life.weight_at_age[a] = T(2 * age - 1);
      life.maturity_at_age[a] = age >= 4 ? T(1.0) : T(0.0);
    }
  }
};

} // namespace bigeye_v2
