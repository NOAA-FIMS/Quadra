#include "../examples/big/catch_at_age_shared.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

int main() {
  example::CatchAtAgeLaplaceModel model;
  for (const auto &row : model.data.age_comp_counts) {
    if (std::accumulate(row.begin(), row.end(), 0) != 200) {
      throw std::runtime_error("composition count row does not sum to 200");
    }
  }

  auto parameters = example::make_big_laplace_parameter_vector();
  const auto fixed_indices = quadra::build_fixed_index(parameters);
  const auto random_indices = quadra::build_random_index(parameters);
  if (fixed_indices.size() != 10 || random_indices.size() != 30) {
    throw std::runtime_error("unexpected big-model parameter partition");
  }

  std::vector<double> full;
  full.reserve(parameters.size());
  for (const auto &parameter : parameters.params) {
    full.push_back(parameter.value);
  }
  const double baseline = model(full);
  full[9] = std::log(10.0);
  const double changed_concentration = model(full);
  if (!std::isfinite(baseline) || !std::isfinite(changed_concentration)) {
    throw std::runtime_error("robust composition objective is not finite");
  }
  if (std::abs(baseline - changed_concentration) < 1.0e-6) {
    throw std::runtime_error("composition concentration does not affect objective");
  }

  std::cout << "PASS: big catch-at-age robust Dirichlet-multinomial\n";
  return 0;
}
