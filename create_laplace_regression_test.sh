#!/usr/bin/env bash
set -euo pipefail

mkdir -p tests/laplace

cat > tests/laplace/test_laplace_gradient_regression.cpp <<'CPP'
#include "../../core/optimizer.hpp"
#include <cmath>
#include <iostream>

using namespace quadra;

struct GaussianREModel {
  AD operator()(std::vector<AD>& p) {
    // theta = p[0], u = p[1]
    const AD& theta = p[0];
    const AD& u = p[1];

    // Joint objective:
    // f(theta,u) = 0.5*(theta - 2)^2 + 0.5*exp(theta)*(u - theta)^2
    //
    // u*(theta) = theta
    // H_uu = exp(theta)
    // Laplace objective = 0.5*(theta - 2)^2 + 0.5*theta - constant
    // gradient = (theta - 2) + 0.5
    return 0.5 * (theta - 2.0) * (theta - 2.0) +
           0.5 * exp(theta) * (u - theta) * (u - theta);
  }
};

int main() {
  ParameterVector params;
  params.add("theta", 1.25, false); // fixed
  params.add("u", 0.0, true);       // random

  GaussianREModel model;

  LaplaceOptions opts = default_laplace_options();
  OptResult result = optimize_lbfgs(model, params, opts);

  const double theta = result.par.at(0);
  const double analytic_expected = (theta - 2.0) + 0.5;
  const double got = result.fixed_gradient.at(0);

  const double diff = std::abs(got - analytic_expected);

  std::cout << "theta," << theta << "\n";
  std::cout << "gradient_got," << got << "\n";
  std::cout << "gradient_expected," << analytic_expected << "\n";
  std::cout << "abs_diff," << diff << "\n";

  if (!(diff < 1e-6)) {
    std::cerr << "FAILED: Laplace profiled gradient regression\n";
    return 1;
  }

  std::cout << "PASSED: Laplace profiled gradient regression\n";
  return 0;
}
CPP

cat > run_laplace_regression_test.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

c++ ${CXXFLAGS:-} -std=c++17 -O2 \
  -I. -Iexternal/eigen \
  tests/laplace/test_laplace_gradient_regression.cpp \
  -o build/tests/test_laplace_gradient_regression

./build/tests/test_laplace_gradient_regression
SH

chmod +x run_laplace_regression_test.sh

echo "created:"
echo "  tests/laplace/test_laplace_gradient_regression.cpp"
echo "  run_laplace_regression_test.sh"
echo
echo "run:"
echo "  ./run_laplace_regression_test.sh"
