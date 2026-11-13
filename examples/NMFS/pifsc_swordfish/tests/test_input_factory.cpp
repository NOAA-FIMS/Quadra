#include "swordfish_input.hpp"
#include <iostream>
/// Run test input factory; return zero on success and nonzero on a failed check
/// or invalid input.
int main() {
  const auto d = imported_swordfish_data();
  size_t included = 0;
  for (const auto &c : d.compositions)
    included += c.use;
  if (d.catches.size() != 2223 || d.indices.size() != 164 ||
      d.compositions.size() != 446 || included != 351 ||
      d.parameters.size() != 74 || d.index_specs.at(24).at(1) != 33 ||
      d.biology.growth[1].length1 != 83.2 ||
      d.compositions.front().sample_size != 17.6 ||
      d.compositions.front().values.at(15) != 1)
    return 1;
  std::cout << "Generated C++ original-observation factory: PASS\n";
}
