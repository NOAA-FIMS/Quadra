#include "../architecture/steps/fleet/logistic_selectivity.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {

bool nearly_equal(double a, double b, double tol = 1.0e-12) {
  return std::abs(a - b) <= tol;
}

} // namespace

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;

  FleetParameters<double> p;
  p.sel_a50 = 5.0;
  p.sel_slope = 1.0;

  FleetState<double> fleet;

  LogisticSelectivity{}(data, p, fleet);

  constexpr double expected[kAges] = {0.01798620996209156, 0.04742587317756678,
                                      0.11920292202211755, 0.26894142136999510,
                                      0.50000000000000000, 0.73105857863000490,
                                      0.88079707797788230, 0.95257412682243340,
                                      0.98201379003790850, 0.99330714907571530};

  for (int a = 0; a < kAges; ++a) {
    if (!nearly_equal(fleet.selectivity_at_age[a], expected[a])) {
      std::cerr << std::setprecision(17) << "FAIL: CAA selectivity_at_age[" << a
                << "] got " << fleet.selectivity_at_age[a] << " expected "
                << expected[a] << " diff "
                << (fleet.selectivity_at_age[a] - expected[a]) << "\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 CAA logistic selectivity regression\n";
  return EXIT_SUCCESS;
}
