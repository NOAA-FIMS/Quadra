#pragma once

#include "../diagnostics/pollock_fixed_effect_diagnostics.hpp"

#include "../../../../../core/optimizer.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace pollock_example {

inline void write_recruitment_deviations(const std::string &path,
                                         const quadra::OptResult &fit,
                                         double rec_rho_report = 0.60) {
  std::ofstream rec(path);
  rec << "year,log_rec_dev,ar1_rho,innovation\n";

  for (std::size_t i = 0; i < fit.u_hat.size(); ++i) {
    const double innovation =
        (i == 0) ? fit.u_hat[i]
                 : (fit.u_hat[i] - rec_rho_report * fit.u_hat[i - 1]);
    rec << (i + 1) << "," << fit.u_hat[i] << "," << rec_rho_report << ","
        << innovation << "\n";
  }
}

inline void print_fit_and_structure_diagnostics(const quadra::OptResult &fit) {
  std::cout << "\nFit diagnostics\n";
  std::cout << "---------------\n";
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "objective          " << fit.value << "\n";
  std::cout << "grad_norm          " << fit.grad_norm << "\n";
  std::cout << "iterations         " << fit.iterations << "\n";
  std::cout << "converged          " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message            " << fit.message << "\n";

  if (!fit.fixed_gradient.empty()) {
    const std::size_t max_grad_i = max_fixed_gradient_index(fit);
    const std::string max_grad_name =
        (max_grad_i < fit.fixed_gradient_names.size())
            ? fit.fixed_gradient_names[max_grad_i]
            : ("fixed_" + std::to_string(max_grad_i));
    std::cout << "max_grad_param     " << max_grad_name << "\n";
    std::cout << "max_grad_value     " << fit.fixed_gradient[max_grad_i]
              << "\n";
    std::cout << "max_abs_grad       "
              << std::abs(fit.fixed_gradient[max_grad_i]) << "\n";
  }

  std::cout << "\nOptimizer structure diagnostics\n";
  std::cout << "-------------------------------\n";
  std::cout << "random effects     " << fit.pattern.random_effect_count << "\n";
  std::cout << "pattern available  " << (fit.pattern.available ? "yes" : "no")
            << "\n";
  std::cout << "detected structure " << fit.pattern.detected_structure << "\n";
  std::cout << "Hessian nonzeros   " << fit.pattern.nonzeros << "\n";
}

inline void print_output_manifest() {
  std::cout << "\nWrote outputs:\n";
  std::cout << "  "
               "examples/NMFS/afsc_walleye_pollock/outputs/"
               "walleye_pollock_fit_summary.csv\n";
  std::cout << "  "
               "examples/NMFS/afsc_walleye_pollock/outputs/"
               "walleye_pollock_recruitment_deviations.csv\n";
}

} // namespace pollock_example
