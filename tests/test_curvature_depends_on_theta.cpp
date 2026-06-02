
#include "test_common.hpp"

DECLARE_ADGRAPH();

struct CurvatureDependsOnTheta {
  int n_random;

  template <typename T> T operator()(const std::vector<T> &p) const {
    T theta = p[0];

    // Weak prior/penalty keeps the optimum finite:
    // Laplace profile objective should have optimum theta = -0.5 * n_random.
    T nll = 0.5 * theta * theta;

    for (int i = 0; i < n_random; ++i) {
      T u = p[1 + i];
      nll += 0.5 * exp(theta) * u * u;
    }

    return nll;
  }
};

int main() {
  using namespace quadra;

  for (int n_random : {1, 2, 5, 10, 25}) {
    quadra_tests::print_banner("Testing CurvatureDependsOnTheta n_random = " +
                               std::to_string(n_random));

    ParameterVector params;
    params.add({"theta", 0.0, ParameterTransform::Identity, false});

    for (int i = 0; i < n_random; ++i) {
      params.add(
          {"u_" + std::to_string(i), 0.0, ParameterTransform::Identity, true});
    }

    CurvatureDependsOnTheta model{n_random};

    auto opts = quadra_tests::default_test_options();
    auto fit = optimize_lbfgs(model, params, opts);

    double theta_hat = fit.par[0];
    double expected = -0.5 * static_cast<double>(n_random);

    std::cout << "theta_hat = " << theta_hat << "\n";
    std::cout << "expected  = " << expected << "\n";
    std::cout << "fx        = " << fit.value << "\n";

    if (std::abs(theta_hat - expected) > 1e-5) {
      std::cerr << "FAIL: theta_hat differs from expected optimum.\n";
      return 1;
    }
  }

  std::cout << "\nPASS\n";
  return 0;
}
