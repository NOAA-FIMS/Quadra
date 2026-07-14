#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LognormalIndexLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.index_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_index_by_year.size(); ++y) {
      const T obs = data.observed_index_by_year[y];
      const T pred = fleet.predicted_index_by_year[y];
      const T r = (std::log(obs) - std::log(pred)) / parameters.index_sigma;

      likelihood.index_nll +=
          T(0.5) * r * r + std::log(parameters.index_sigma) + std::log(obs);
    }

    likelihood.total_nll += likelihood.index_nll;
  }
};

} // namespace bigeye_v2
