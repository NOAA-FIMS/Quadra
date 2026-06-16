#pragma once

#include <string>
#include <vector>

namespace pollock_example {

inline std::string pollock_output_dir() {
  return "examples/NMFS/afsc_walleye_pollock/outputs";
}

inline std::string pollock_output_path(const std::string &filename) {
  return pollock_output_dir() + "/" + filename;
}

inline std::vector<std::string> pollock_random_effect_names(std::size_t n) {
  std::vector<std::string> names;
  names.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    names.push_back("log_rec_dev_" + std::to_string(i + 1));
  }
  return names;
}

} // namespace pollock_example
