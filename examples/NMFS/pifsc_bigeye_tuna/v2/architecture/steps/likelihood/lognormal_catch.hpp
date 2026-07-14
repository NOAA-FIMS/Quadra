#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LognormalCatchLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.catch_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_catch_biomass_by_year.size();
         ++y) {
      const T obs = data.observed_catch_biomass_by_year[y];
      const T pred = fleet.total_catch_biomass_by_year[y];
      const T r = (std::log(obs) - std::log(pred)) / parameters.catch_sigma;

      likelihood.catch_nll +=
          T(0.5) * r * r + std::log(parameters.catch_sigma) + std::log(obs);
    }

    likelihood.total_nll += likelihood.catch_nll;
  }
};

} // namespace bigeye_v2
