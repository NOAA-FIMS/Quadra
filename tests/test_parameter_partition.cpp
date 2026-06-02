#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../core/model/parameter_partition.hpp"

int main() {
  std::cout << "Testing fixed/random parameter partitioning\n";

  quadra::ParameterSet parameters;
  parameters.add("log_R0", 10.0, quadra::ParameterTransform::Identity, false);
  parameters.add("log_sigma_R", -1.0, quadra::ParameterTransform::Log, false);
  parameters.add("u_1", 0.1, quadra::ParameterTransform::Identity, true);
  parameters.add("u_2", -0.2, quadra::ParameterTransform::Identity, true);
  parameters.add("log_q", -3.0, quadra::ParameterTransform::Log, false);

  const auto partition = quadra::partition_parameters(parameters);

  if (partition.n_fixed() != 3) {
    std::cerr << "FAIL: expected 3 fixed effects\n";
    return 1;
  }

  if (partition.n_random() != 2) {
    std::cerr << "FAIL: expected 2 random effects\n";
    return 1;
  }

  const std::vector<size_t> expected_fixed = {0, 1, 4};
  const std::vector<size_t> expected_random = {2, 3};

  if (partition.fixed_indices_m != expected_fixed) {
    std::cerr << "FAIL: fixed indices mismatch\n";
    return 1;
  }

  if (partition.random_indices_m != expected_random) {
    std::cerr << "FAIL: random indices mismatch\n";
    return 1;
  }

  const auto fixed_names = quadra::fixed_effect_names(parameters);
  const auto random_names = quadra::random_effect_names(parameters);

  if (fixed_names !=
      std::vector<std::string>{"log_R0", "log_sigma_R", "log_q"}) {
    std::cerr << "FAIL: fixed names mismatch\n";
    return 1;
  }

  if (random_names != std::vector<std::string>{"u_1", "u_2"}) {
    std::cerr << "FAIL: random names mismatch\n";
    return 1;
  }

  std::vector<double> full = {10.0, -1.0, 0.1, -0.2, -3.0};
  auto split = quadra::split_parameters(full, partition);

  if (split.fixed_m != std::vector<double>{10.0, -1.0, -3.0}) {
    std::cerr << "FAIL: fixed split mismatch\n";
    return 1;
  }

  if (split.random_m != std::vector<double>{0.1, -0.2}) {
    std::cerr << "FAIL: random split mismatch\n";
    return 1;
  }

  auto merged = quadra::merge_parameters(split, partition);

  if (merged != full) {
    std::cerr << "FAIL: merged vector did not round-trip\n";
    return 1;
  }

  std::cout << "fixed indices:";
  for (size_t i : partition.fixed_indices_m) {
    std::cout << " " << i;
  }
  std::cout << "\n";

  std::cout << "random indices:";
  for (size_t i : partition.random_indices_m) {
    std::cout << " " << i;
  }
  std::cout << "\n";

  std::cout << "PASS\n";
  return 0;
}
