#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"
#include "../../steps/likelihood/lognormal_catch.hpp"
#include "../../steps/likelihood/lognormal_index.hpp"
#include "../../steps/likelihood/multinomial_agecomp.hpp"
#include "likelihood_context.hpp"

namespace bigeye_v2 {

struct LikelihoodPackage {

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const LikelihoodContext<T> &context) const {
    (*this)(data, *context.parameters, *context.fleet, *context.likelihood);
  }

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.catch_nll = T(0.0);
    likelihood.index_nll = T(0.0);
    likelihood.agecomp_nll = T(0.0);
    likelihood.total_nll = T(0.0);

    LognormalCatchLikelihood{}(data, parameters, fleet, likelihood);
    LognormalIndexLikelihood{}(data, parameters, fleet, likelihood);
    MultinomialAgeCompLikelihood{}(data, fleet, likelihood);
  }
};

} // namespace bigeye_v2
