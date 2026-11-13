#include "../core/had_quadra.hpp"
#include <iostream>
#include <stdexcept>

/// Verify value-based branches independently of directional derivative seeds.
/// Equal primal values must compare equal even when their seeds differ; a
/// comparison chooses a branch but must not mutate either operand's derivatives.
int main() {
  const had::ThirdOrderScalar low(1.0, 7.0), high(2.0, -3.0), same(1.0, -9.0);
  if (!(low < high) || !(low <= high) || !(high > low) || !(high >= low) ||
      !(low != high) || low == high || !(low == same) || low != same ||
      !(low <= same) || !(low >= same) || low < same || low > same ||
      low.d1 != 7.0 || high.d1 != -3.0 || same.d1 != -9.0)
    throw std::runtime_error("ThirdOrderScalar comparison semantics changed");
  std::cout << "Third-order scalar value comparisons: PASS\n";
}
