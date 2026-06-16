#pragma once

#include <algorithm>
#include <cstddef>
#include <string>

#include "../../../../../core/optimizer.hpp"

namespace pollock {

inline quadra::ParameterVector make_params(std::size_t n_years) {
  quadra::ParameterVector p;

  auto add_param = [&](const std::string &name, double value, bool random) {
    p.add(quadra::Parameter(
        name,
        value,
        quadra::ParameterTransform::Identity,
        random));
  };

  add_param("log_r0", 8.0, false);
  add_param("log_fbar", -3.7, false);

#ifdef WALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT
  const std::size_t n_random_recruitment =
      std::min<std::size_t>(
          n_years,
          static_cast<std::size_t>(WALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT));

  for (std::size_t i = 0; i < n_random_recruitment; ++i) {
    add_param("log_rec_dev_" + std::to_string(i + 1), 0.0, true);
  }
#else
  (void)n_years;
#endif

  return p;
}

}  // namespace pollock
