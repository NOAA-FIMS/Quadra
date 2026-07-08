#include "../architecture/packages/movement/movement_package.hpp"

#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 2;

  MovementParameters<double> parameters;

  std::vector<PopulationState<double>> populations(2);
  for (auto &pop : populations) {
    pop.numbers_at_age.assign(
        static_cast<std::size_t>(data.n_years),
        std::array<double, kAges>{});
  }

  populations[0].numbers_at_age[0][0] = 1000.0;
  populations[1].numbers_at_age[0][0] = 500.0;

  MovementContext<double> context{&parameters, &populations};

  MovementPackage{}(data, context);

  if (populations[0].numbers_at_age[0][0] != 1000.0 ||
      populations[1].numbers_at_age[0][0] != 500.0) {
    std::cerr << "FAIL: identity movement changed population numbers\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA MovementPackage regression\n";
  return EXIT_SUCCESS;
}
