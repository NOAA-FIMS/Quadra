
#include <cmath>
#include <iostream>
#include <vector>

#include "../core/optimizer.hpp"

DECLARE_ADGRAPH();

// CPUE index model with latent abundance deviations.
// This is a minimal fisheries-relevant Laplace example, not a full assessment.

struct CPUEWithLatentAbundance {
  std::vector<double> log_cpue;

  template <typename T> T operator()(const std::vector<T> &p) const {
    T log_q = p[0];
    T log_sigma_obs = p[1];
    T log_sigma_proc = p[2];

    T sigma_obs = exp(log_sigma_obs);
    T sigma_proc = exp(log_sigma_proc);

    T nll = 0.0;

    for (int t = 0; t < static_cast<int>(log_cpue.size()); ++t) {
      T logN_t = p[3 + t];

      T pred = log_q + logN_t;
      T r_obs = (log_cpue[t] - pred) / sigma_obs;
      nll += 0.5 * r_obs * r_obs + log_sigma_obs;

      if (t == 0) {
        T z = logN_t / sigma_proc;
        nll += 0.5 * z * z + log_sigma_proc;
      } else {
        T logN_prev = p[3 + t - 1];
        T z = (logN_t - logN_prev) / sigma_proc;
        nll += 0.5 * z * z + log_sigma_proc;
      }
    }

    nll += 0.001 * log_q * log_q;
    nll += 0.001 * log_sigma_obs * log_sigma_obs;
    nll += 0.001 * log_sigma_proc * log_sigma_proc;

    return nll;
  }
};

int main() {
  using namespace quadra;

  std::vector<double> log_cpue = {
      std::log(1.0),  std::log(0.95), std::log(0.88), std::log(0.80),
      std::log(0.77), std::log(0.70), std::log(0.66), std::log(0.62)};

  ParameterVector params;
  params.add({"log_q", 0.0, ParameterTransform::Identity, false});
  params.add(
      {"log_sigma_obs", std::log(0.2), ParameterTransform::Identity, false});
  params.add(
      {"log_sigma_proc", std::log(0.15), ParameterTransform::Identity, false});

  for (int t = 0; t < static_cast<int>(log_cpue.size()); ++t) {
    params.add(
        {"logN_" + std::to_string(t), 0.0, ParameterTransform::Identity, true});
  }

  CPUEWithLatentAbundance model{log_cpue};

  LaplaceOptions opts;
  opts.use_hutchinson_trace = false;
  opts.hessian_drop_tol = 0.0;

  auto fit = optimize_lbfgs(model, params, opts);

  std::cout << "fit.value = " << fit.value << "\n";
  std::cout << "log_q_hat = " << fit.par[0] << "\n";

  return 0;
}
