#pragma once

#include "../../../../../core/optimizer.hpp"

#include <fstream>
#include <iomanip>
#include <string>

namespace pollock_example {

void write_summary(const std::string &path, const quadra::OptResult &fit)
{
  std::ofstream out(path);
  out << std::setprecision(15);
  out << "field,value\n";
  out << "objective," << fit.value << "\n";
  out << "grad_norm," << fit.grad_norm << "\n";
  out << "iterations," << fit.iterations << "\n";
  out << "converged," << (fit.converged ? "yes" : "no") << "\n";
  out << "message," << fit.message << "\n";
  out << "random_effects," << fit.u_hat.size() << "\n";
}


}  // namespace pollock_example

// Compatibility alias for current walleye_pollock.cpp call sites.
using pollock_example::write_summary;
