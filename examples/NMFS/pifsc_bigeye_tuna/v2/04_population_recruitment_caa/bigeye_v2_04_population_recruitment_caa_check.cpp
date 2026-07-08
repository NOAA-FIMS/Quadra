#include "../architecture/parameters/population_parameters.hpp"
#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/population/recruitment.hpp"
#include "../common/model_data.hpp"

#include <cstdlib>
#include <iostream>

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 4;

  PopulationParameters<double> p;
  p.r0 = 1234.5;

  PopulationState<double> population;

  FixedRecruitment{}(data, p, population);

  if (population.recruits_by_year.size() != 4) {
    std::cerr << "FAIL: recruits size\n";
    return EXIT_FAILURE;
  }

  for (double r : population.recruits_by_year) {
    if (r != 1234.5) {
      std::cerr << "FAIL: recruit value\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 CAA population recruitment regression\n";
  return EXIT_SUCCESS;
}
