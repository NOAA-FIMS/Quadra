#pragma once

#include "../../../common/bigeye_constants.hpp"
#include "../../../common/model_data.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"

#include <cmath>

namespace bigeye_v2 {

struct MultinomialAgeCompLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data, const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.agecomp_nll = T(0.0);

    const T eps = T(1.0e-12);

    for (std::size_t y = 0; y < data.observed_catch_age_proportion.size();
         ++y) {
      const T n_eff = data.catch_agecomp_sample_size[y];

      for (int a = 0; a < kAges; ++a) {
        const T obs = data.observed_catch_age_proportion[y][a];
        const T pred = fleet.predicted_catch_age_proportion[y][a] + eps;

        likelihood.agecomp_nll -= n_eff * obs * std::log(pred);
      }
    }

    likelihood.total_nll += likelihood.agecomp_nll;
  }
};

} // namespace bigeye_v2
