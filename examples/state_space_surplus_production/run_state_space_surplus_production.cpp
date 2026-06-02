#include "state_space_surplus_production.hpp"

#include <cmath>
#include <iostream>

namespace ss = quadra_examples::state_space_surplus_production;

int main() {
  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();
  const std::vector<double> u = ss::zero_random_effects(data);

  ss::print_report(data, par, u);

  return 0;
}
