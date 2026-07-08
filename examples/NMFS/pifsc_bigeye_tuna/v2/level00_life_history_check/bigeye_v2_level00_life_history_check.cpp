#include "../common/life_history.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

bool nearly_equal(double a, double b, double tol = 1.0e-12) {
  return std::abs(a - b) <= tol;
}

} // namespace

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  BigeyeModelParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);

  BigeyeDerived<double> d;

  BigeyeLifeHistory life_history;
  life_history(data, p, d);

  constexpr double expected_m[kAges] = {
      0.3375, 0.3375, 0.3375, 0.45, 0.45,
      0.45,   0.45,   0.2925, 0.2925, 0.2925};

  constexpr double expected_weight[kAges] = {
      1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};

  constexpr double expected_maturity[kAges] = {
      0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  for (int a = 0; a < kAges; ++a) {
    if (!nearly_equal(d.m_at_age[a], expected_m[a])) {
      std::cerr << "FAIL: m_at_age[" << a << "]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.weight_at_age[a], expected_weight[a])) {
      std::cerr << "FAIL: weight_at_age[" << a << "]\n";
      return EXIT_FAILURE;
    }

    if (!nearly_equal(d.maturity_at_age[a], expected_maturity[a])) {
      std::cerr << "FAIL: maturity_at_age[" << a << "]\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 Level00 life-history regression\n";
  return EXIT_SUCCESS;
}
