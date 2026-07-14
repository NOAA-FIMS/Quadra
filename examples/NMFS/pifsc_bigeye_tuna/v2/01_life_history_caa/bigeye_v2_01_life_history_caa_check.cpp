#include "../architecture/steps/life_history/life_history.hpp"

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

  LifeHistoryParameters<double> p;
  p.log_m_young_offset = std::log(0.75);
  p.log_m_old_offset = std::log(0.65);

  LifeHistoryState<double> life;

  BigeyeLifeHistory{}(data, p, life);

  constexpr double expected_m_age_1 = 0.3375;
  constexpr double expected_m_age_4 = 0.45;
  constexpr double expected_m_age_8 = 0.2925;

  if (!nearly_equal(life.m_at_age[0], expected_m_age_1)) {
    std::cerr << std::setprecision(17) << "FAIL: m age 1 got "
              << life.m_at_age[0] << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.m_at_age[3], expected_m_age_4)) {
    std::cerr << std::setprecision(17) << "FAIL: m age 4 got "
              << life.m_at_age[3] << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.m_at_age[7], expected_m_age_8)) {
    std::cerr << std::setprecision(17) << "FAIL: m age 8 got "
              << life.m_at_age[7] << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.weight_at_age[9], 19.0)) {
    std::cerr << "FAIL: weight age 10\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(life.maturity_at_age[2], 0.0) ||
      !nearly_equal(life.maturity_at_age[3], 1.0)) {
    std::cerr << "FAIL: maturity knife-edge\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA life-history regression\n";
  return EXIT_SUCCESS;
}
