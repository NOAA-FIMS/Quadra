#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

namespace bigeye_v2 {

struct FishingMortality {
  template <typename T>
  void operator()(const BigeyeModelData<T> & /*data*/,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    for (int a = 0; a < kAges; ++a) {
      d.f_at_age[a] = p.fbar * d.selectivity_at_age[a];
      d.z_at_age[a] = d.m_at_age[a] + d.f_at_age[a];
    }
  }
};

} // namespace bigeye_v2
