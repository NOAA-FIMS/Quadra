#include <iostream>
#include <vector>

#include "../core/model/parameter.hpp"

int main() {
  quadra::ParameterSet parameters;
  parameters.add("log_sigma", 0.25, quadra::ParameterTransform::Log);
  parameters.add("logit_p", 0.0, quadra::ParameterTransform::Logit);
  parameters.add("sqrt_q", 2.0, quadra::ParameterTransform::Square);

  const auto x = parameters.initials();
  const auto transforms = parameters.transforms();

  const auto constrained = quadra::apply_transforms(x, transforms);
  const auto log_j_terms = quadra::transform_log_jacobians(x, transforms);
  const double total_log_j = quadra::sum_transform_log_jacobian(x, transforms);

  const auto names = parameters.names();

  std::cout << "transform Jacobian demo\n";

  for (size_t i = 0; i < x.size(); ++i) {
    std::cout << names[i] << ": unconstrained = " << x[i]
              << ", constrained = " << constrained[i]
              << ", logJ = " << log_j_terms[i] << "\n";
  }

  std::cout << "total logJ = " << total_log_j << "\n";

  return 0;
}
