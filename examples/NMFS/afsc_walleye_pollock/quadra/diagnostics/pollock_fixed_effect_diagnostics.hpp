#pragma once

#include "../../../../../core/optimizer.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <string>

namespace pollock_example {

void write_fixed_gradient_diagnostics(const std::string &path,
                                    const quadra::OptResult &fit) {
std::ofstream out(path);
out << std::setprecision(15);
out << "parameter,gradient,abs_gradient\n";

for (std::size_t i = 0; i < fit.fixed_gradient.size(); ++i) {
  const std::string name =
      (i < fit.fixed_gradient_names.size()) ? fit.fixed_gradient_names[i]
                                            : ("fixed_" + std::to_string(i));
  const double g = fit.fixed_gradient[i];
  out << name << "," << g << "," << std::abs(g) << "\n";
}
}

std::size_t max_fixed_gradient_index(const quadra::OptResult &fit) {
std::size_t best = 0;
double best_abs = -1.0;

for (std::size_t i = 0; i < fit.fixed_gradient.size(); ++i) {
  const double a = std::abs(fit.fixed_gradient[i]);
  if (a > best_abs) {
    best = i;
    best_abs = a;
  }
}

return best;
}

void write_fixed_parameter_estimates(const std::string &path,
                                   const quadra::OptResult &fit) {
std::ofstream out(path);
out << std::setprecision(15);
out << "parameter,estimate,exp_estimate\n";

for (std::size_t i = 0; i < fit.par.size(); ++i) {
  const std::string name =
      (i < fit.fixed_gradient_names.size()) ? fit.fixed_gradient_names[i]
                                            : ("fixed_" + std::to_string(i));
  out << name << "," << fit.par[i] << "," << std::exp(fit.par[i]) << "\n";
}
}

}  // namespace pollock_example

// Compatibility aliases for the current Pollock implementation, which still
// calls these helpers unqualified from walleye_pollock.cpp.
using pollock_example::max_fixed_gradient_index;
using pollock_example::write_fixed_gradient_diagnostics;
using pollock_example::write_fixed_parameter_estimates;
