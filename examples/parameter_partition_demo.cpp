#include <iostream>
#include <vector>

#include "../core/model/parameter_partition.hpp"

int main() {
  quadra::ParameterSet parameters;
  parameters.add("log_R0", 10.0, quadra::ParameterTransform::Identity, false);
  parameters.add("log_sigma_R", -1.0, quadra::ParameterTransform::Log, false);
  parameters.add("u_1", 0.1, quadra::ParameterTransform::Identity, true);
  parameters.add("u_2", -0.2, quadra::ParameterTransform::Identity, true);
  parameters.add("log_q", -3.0, quadra::ParameterTransform::Log, false);

  const auto partition = quadra::partition_parameters(parameters);
  const auto full = parameters.initials();
  const auto split = quadra::split_parameters(full, partition);

  const auto fixed_names = quadra::fixed_effect_names(parameters);
  const auto random_names = quadra::random_effect_names(parameters);

  std::cout << "fixed effects:\n";
  for (size_t i = 0; i < split.fixed_m.size(); ++i) {
    std::cout << "  " << fixed_names[i] << " = " << split.fixed_m[i] << "\n";
  }

  std::cout << "random effects:\n";
  for (size_t i = 0; i < split.random_m.size(); ++i) {
    std::cout << "  " << random_names[i] << " = " << split.random_m[i] << "\n";
  }

  const auto merged = quadra::merge_parameters(split, partition);

  std::cout << "merged full vector:\n";
  for (double x : merged) {
    std::cout << "  " << x << "\n";
  }

  return 0;
}
