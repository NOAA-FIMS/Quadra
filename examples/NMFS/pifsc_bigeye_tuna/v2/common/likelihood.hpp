#pragma once

#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LognormalCatchLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.catch_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_catch_biomass_by_year.size(); ++y) {
      const T obs = data.observed_catch_biomass_by_year[y];
      const T pred = d.total_catch_biomass_by_year[y];

      const T r = (std::log(obs) - std::log(pred)) / p.catch_sigma;

      d.catch_nll +=
          T(0.5) * r * r +
          std::log(p.catch_sigma) +
          std::log(obs);
    }

    d.total_nll += d.catch_nll;
  }
};


struct LognormalIndexLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.index_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_index_by_year.size(); ++y) {
      const T obs = data.observed_index_by_year[y];
      const T pred = d.predicted_index_by_year[y];

      const T r = (std::log(obs) - std::log(pred)) / p.index_sigma;

      d.index_nll +=
          T(0.5) * r * r +
          std::log(p.index_sigma) +
          std::log(obs);
    }

    d.total_nll += d.index_nll;
  }
};


struct MultinomialAgeCompLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &,
                  BigeyeDerived<T> &d) const {
    d.agecomp_nll = T(0.0);

    const T eps = T(1.0e-12);

    for (std::size_t y = 0; y < data.observed_catch_age_proportion.size(); ++y) {
      const T n_eff = data.catch_agecomp_sample_size[y];

      for (int a = 0; a < kAges; ++a) {
        const T obs = data.observed_catch_age_proportion[y][a];
        const T pred = d.predicted_catch_age_proportion[y][a] + eps;

        d.agecomp_nll -= n_eff * obs * std::log(pred);
      }
    }

    d.total_nll += d.agecomp_nll;
  }
};

} // namespace bigeye_v2
