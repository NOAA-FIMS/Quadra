#include "../include/quadra/quadra.hpp"

#include <iostream>

int main()
{
    quadra::LaplaceImplicitWorkspace workspace;
    quadra::LaplaceProfiledDerivedReport report;
    quadra::FixedEffectCovarianceResult covariance;

    (void)workspace;
    (void)report;
    (void)covariance;

    std::cout << "PASS: public Quadra API headers compile\n";

    return 0;
}
