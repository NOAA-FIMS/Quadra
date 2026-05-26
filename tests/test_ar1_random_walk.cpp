
#include "test_common.hpp"

DECLARE_ADGRAPH();

struct GaussianRandomWalk {
  std::vector<double> y;

  template <typename T> T operator()(const std::vector<T> &p) const {
    T mu = p[0];
    T log_sigma_obs = p[1];
    T log_sigma_rw = p[2];

    T sigma_obs = exp(log_sigma_obs);
    T sigma_rw = exp(log_sigma_rw);

    T nll = 0.0;

    for (int t = 0; t < static_cast<int>(y.size()); ++t) {
      T x_t = p[3 + t];

      // Observation model
      T r = (y[t] - (mu + x_t)) / sigma_obs;
      nll += 0.5 * r * r + log_sigma_obs;

      // Random walk prior
      if (t == 0) {
        T z = x_t / sigma_rw;
        nll += 0.5 * z * z + log_sigma_rw;
      } else {
        T x_prev = p[3 + t - 1];
        T z = (x_t - x_prev) / sigma_rw;
        nll += 0.5 * z * z + log_sigma_rw;
      }
    }

    // weak fixed-effect penalties
    nll += 0.001 * mu * mu;
    nll += 0.001 * log_sigma_obs * log_sigma_obs;
    nll += 0.001 * log_sigma_rw * log_sigma_rw;

    return nll;
  }
};

int main() {
  using namespace quadra;
  quadra_tests::print_banner("Testing GaussianRandomWalk");

  std::vector<double> y = {1.0, 1.2, 0.9, 1.5, 1.8, 1.7, 2.0, 2.2};

  ParameterVector params;
  params.add({"mu", 1.0, ParameterTransform::Identity, false});
  params.add(
      {"log_sigma_obs", std::log(0.4), ParameterTransform::Identity, false});
  params.add(
      {"log_sigma_rw", std::log(0.3), ParameterTransform::Identity, false});

  for (int t = 0; t < static_cast<int>(y.size()); ++t) {
    params.add(
        {"x_" + std::to_string(t), 0.0, ParameterTransform::Identity, true});
  }

  GaussianRandomWalk model{y};

  auto opts = quadra_tests::default_test_options();
  auto fit = optimize_lbfgs(model, params, opts);

  std::cout << "fit.value = " << fit.value << "\n";
  std::cout << "mu_hat = " << fit.par[0] << "\n";

  if (!std::isfinite(fit.value)) {
    std::cerr << "FAIL: non-finite objective.\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
