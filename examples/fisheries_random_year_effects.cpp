
#include <cmath>
#include <iostream>
#include <vector>

#include "../core/optimizer.hpp"

DECLARE_ADGRAPH();

// Simple fisheries-style biomass index model with year random effects.
//
// log I_t ~ Normal(log(q) + log(B_t) + u_t, sigma_obs)
// u_t follows iid normal prior.
// This is not a full stock assessment model; it is a compact Laplace example.

struct FisheriesRandomYearEffects {
  std::vector<double> log_index;
  std::vector<double> log_biomass;

  template <typename T> T operator()(const std::vector<T> &p) const {
    T log_q = p[0];
    T log_sigma_obs = p[1];
    T log_sigma_re = p[2];

    T sigma_obs = exp(log_sigma_obs);
    T sigma_re = exp(log_sigma_re);

    T nll = 0.0;

    for (int t = 0; t < static_cast<int>(log_index.size()); ++t) {
      T u_t = p[3 + t];

      T pred = log_q + log_biomass[t] + u_t;
      T r = (log_index[t] - pred) / sigma_obs;

      nll += 0.5 * r * r + log_sigma_obs;

      T z = u_t / sigma_re;
      nll += 0.5 * z * z + log_sigma_re;
    }

    // weak fixed-effect penalties
    nll += 0.001 * log_q * log_q;
    nll += 0.001 * log_sigma_obs * log_sigma_obs;
    nll += 0.001 * log_sigma_re * log_sigma_re;

    return nll;
  }
};

int main() {
  using namespace quadra;

  std::vector<double> log_biomass = {
      std::log(1000), std::log(950), std::log(900), std::log(860),
      std::log(830),  std::log(790), std::log(760), std::log(720)};

  std::vector<double> log_index = {std::log(10.1), std::log(9.4), std::log(9.2),
                                   std::log(8.5),  std::log(8.1), std::log(7.8),
                                   std::log(7.6),  std::log(7.1)};

  ParameterVector params;
  params.add({"log_q", std::log(0.01), ParameterTransform::Identity, false});
  params.add(
      {"log_sigma_obs", std::log(0.2), ParameterTransform::Identity, false});
  params.add(
      {"log_sigma_re", std::log(0.1), ParameterTransform::Identity, false});

  for (int t = 0; t < static_cast<int>(log_index.size()); ++t) {
    params.add({"u_year_" + std::to_string(t), 0.0,
                ParameterTransform::Identity, true});
  }

  FisheriesRandomYearEffects model{log_index, log_biomass};

  LaplaceOptions opts;
  opts.use_hutchinson_trace = false; // small example: exact deterministic trace
  opts.hessian_drop_tol = 0.0;

  auto fit = optimize_lbfgs(model, params, opts);

  // Optimizer result diagnostics.
  //
  // The example does not manually inspect Hessians or select a backend.
  // Quadra returns these diagnostics through the public fit object.
  std::cout << "\nOptimizer diagnostics\n";
  std::cout << "---------------------\n";
  std::cout << "Converged:          " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "Message:            " << fit.message << "\n";
  std::cout << "Random effects:     " << fit.pattern.random_effect_count << "\n";
  std::cout << "Pattern available:  " << (fit.pattern.available ? "yes" : "no") << "\n";
  std::cout << "Detected structure: " << fit.pattern.detected_structure << "\n";
  std::cout << "Laplace backend:    " << fit.pattern.backend << "\n";
  std::cout << "Random solver:      " << fit.pattern.solver << "\n";
  std::cout << "Expected complexity:" << fit.pattern.complexity << "\n";
  std::cout << "Bandwidth:          " << fit.pattern.bandwidth << "\n";
  std::cout << "Hessian nonzeros:   " << fit.pattern.nonzeros << "\n";

  std::cout << "fit.value = " << fit.value << "\n";
  std::cout << "log_q_hat = " << fit.par[0] << "\n";

  return 0;
}
