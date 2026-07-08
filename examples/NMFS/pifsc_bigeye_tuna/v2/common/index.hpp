#pragma once

#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

namespace bigeye_v2 {

struct BiomassIndexPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.predicted_index_by_year.assign(
        static_cast<std::size_t>(data.n_years), T(0.0));

    for (std::size_t y = 0; y < d.spawning_biomass_by_year.size(); ++y) {
      d.predicted_index_by_year[y] = p.q_index * d.spawning_biomass_by_year[y];
    }
  }
};

} // namespace bigeye_v2
