#include "swordfish/benchmarks.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish_input.hpp"
#include <iostream>
/// Run test benchmarks; return zero on success and nonzero on a failed check or
/// invalid input.
int main() {
  auto data = imported_swordfish_data();
  auto p = swordfish::read_reference_parameters(
      data, "build/biology_fitted_parameters.txt",
      "build/population_fitted_deviations.txt");
  swordfish::Benchmarks b(data, p.first, p.second);
  const auto virgin = b(0), spr = b.target(.4, false),
             target = b.target(.3, true), msy = b.msy();
  if (std::abs(virgin.spr - 1) > 1e-12 || virgin.yield != 0 ||
      std::abs(virgin.annual_f) > 1e-12 || std::abs(spr.spr - .4) > 1e-12 ||
      std::abs(target.ssb / virgin.ssb - .3) > 1e-12 ||
      msy.yield < b(msy.scale * .99).yield ||
      msy.yield < b(msy.scale * 1.01).yield)
    throw std::runtime_error("Benchmark root or maximum failure");
  std::cout << "Equilibrium roots and MSY maximum pass\n";
}
