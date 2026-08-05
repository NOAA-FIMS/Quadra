#include "../include/quadra/quadra.hpp"

#include <iostream>

int main() {
  quadra::LaplaceImplicitWorkspace workspace;
  quadra::LaplaceProfiledDerivedReport report;
  quadra::FixedEffectCovarianceResult covariance;
  const double standard_normal_at_zero =
      quadra::stats::normal_logpdf(0.0, 0.0, 1.0);

  (void)workspace;
  (void)report;
  (void)covariance;
  (void)standard_normal_at_zero;

  std::cout << "PASS: public Quadra API headers compile\n";

  return 0;
}
