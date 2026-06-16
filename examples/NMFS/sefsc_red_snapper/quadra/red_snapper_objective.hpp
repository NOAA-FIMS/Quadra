#pragma once

#include "red_snapper_age_structured.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace sefsc_red_snapper {

struct ObjectiveOptions {
  double sigma_log_index = 0.20;
  double sigma_log_catch = 0.15;
  double min_positive = 1.0e-12;
};

struct ObjectiveBreakdown {
  double total = 0.0;
  double index_nll = 0.0;
  double catch_nll = 0.0;
  int n_index = 0;
  int n_catch = 0;
};

inline double square(double x) { return x * x; }

inline double lognormal_nll_no_constant(double observed, double predicted,
                                        double sigma, double min_positive) {
  const double obs = std::max(observed, min_positive);
  const double pred = std::max(predicted, min_positive);
  const double z = (std::log(obs) - std::log(pred)) / sigma;
  return 0.5 * square(z);
}

inline ObjectiveBreakdown evaluate_objective_breakdown(
    const std::vector<Observation> &observations,
    const AgeStructuredParams &params,
    const ObjectiveOptions &options = ObjectiveOptions{}) {
  ObjectiveBreakdown out;

  const auto rows =
      run_deterministic_age_structured_model(observations, params);
  if (rows.size() != observations.size()) {
    throw std::runtime_error("Objective trajectory/observation size mismatch");
  }

  for (std::size_t i = 0; i < observations.size(); ++i) {
    const auto &obs = observations[i];
    const auto &pred = rows[i];

    if (std::isfinite(obs.index) && obs.index > 0.0) {
      const double nll = lognormal_nll_no_constant(obs.index, pred.index_hat,
                                                   options.sigma_log_index,
                                                   options.min_positive);
      out.index_nll += nll;
      ++out.n_index;
    }

    if (std::isfinite(obs.catch_mt) && obs.catch_mt > 0.0) {
      const double nll = lognormal_nll_no_constant(obs.catch_mt, pred.catch_hat,
                                                   options.sigma_log_catch,
                                                   options.min_positive);
      out.catch_nll += nll;
      ++out.n_catch;
    }
  }

  out.total = out.index_nll + out.catch_nll;
  return out;
}

inline double
evaluate_objective(const std::vector<Observation> &observations,
                   const AgeStructuredParams &params,
                   const ObjectiveOptions &options = ObjectiveOptions{}) {
  return evaluate_objective_breakdown(observations, params, options).total;
}

} // namespace sefsc_red_snapper
