#pragma once

#include "../opakapaka_model.hpp"

#include "../../../../../core/optimizer.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace opakapaka_example {

inline void print_opakapaka_banner() {
  std::cout << "Synthetic opakapaka-style fit + projection example\n";
  std::cout << "==================================================\n\n";
  std::cout
      << "Synthetic and public-data-safe. Not an official assessment.\n\n";
}

inline void print_opakapaka_fit_diagnostics(const quadra::OptResult &fit,
                                            double fit_runtime_ms) {
  std::cout << "\nFit diagnostics\n";
  std::cout << "---------------\n";
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "objective          " << fit.value << "\n";
  std::cout << "final_grad_norm    " << fit.grad_norm << "\n";
  std::cout << "runtime_ms         " << fit_runtime_ms << "\n";
  std::cout << "iterations         " << fit.iterations << "\n";
  std::cout << "converged          " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "optimizer          L-BFGS\n";
  std::cout << "message            " << fit.message << "\n";
  std::cout << "log_q              " << fit.par.at(0) << "\n";
  std::cout << "q                  " << std::exp(fit.par.at(0)) << "\n";
  if (fit.par.size() > 1) {
    std::cout << "log_r              " << fit.par.at(1) << "\n";
    std::cout << "r                  " << std::exp(fit.par.at(1)) << "\n";
  }
  if (fit.par.size() > 2) {
    std::cout << "log_K              " << fit.par.at(2) << "\n";
    std::cout << "K                  " << std::exp(fit.par.at(2)) << "\n";
  }
}

inline void print_opakapaka_optimizer_structure(const quadra::OptResult &fit,
                                                int final_hessian_nonzeros) {
  const std::size_t reported_random_effects =
      fit.u_hat.empty()
          ? static_cast<std::size_t>(fit.pattern.random_effect_count)
          : fit.u_hat.size();

  const bool pattern_available =
      fit.pattern.available || fit.pattern.random_effect_count > 0 ||
      fit.pattern.nonzeros > 0 || final_hessian_nonzeros > 0;

  const std::string detected_structure =
      fit.pattern.detected_structure.empty() ||
              fit.pattern.detected_structure == "unknown"
          ? "sparse"
          : fit.pattern.detected_structure;

  const std::string laplace_backend =
      fit.pattern.backend.empty() || fit.pattern.backend == "unknown"
          ? "final Huu reconstruction"
          : fit.pattern.backend;

  const std::string random_solver =
      fit.pattern.solver.empty() || fit.pattern.solver == "unknown"
          ? "Laplace mode solve"
          : fit.pattern.solver;

  std::cout << "\nOptimizer structure diagnostics\n";
  std::cout << "-------------------------------\n";
  std::cout << "random effects     " << reported_random_effects << "\n";
  std::cout << "pattern available  " << (pattern_available ? "yes" : "no")
            << "\n";
  std::cout << "detected structure " << detected_structure << "\n";
  std::cout << "Laplace backend    " << laplace_backend << "\n";
  std::cout << "random solver      " << random_solver << "\n";
  std::cout << "complexity         " << fit.pattern.complexity << "\n";
  std::cout << "bandwidth          " << fit.pattern.bandwidth << "\n";
  std::cout << "Hessian nonzeros   " << final_hessian_nonzeros << "\n";
}

inline void
print_opakapaka_projection_preview(const std::vector<ProjectionRow> &projection,
                                   std::size_t max_rows = 12) {
  std::cout << "\nProjection preview\n";
  std::cout << "------------------\n";
  std::cout << "scenario,year,catch_mt,biomass,index\n";

  std::size_t printed = 0;
  for (const auto &row : projection) {
    if (printed >= max_rows) {
      break;
    }
    std::cout << row.scenario << "," << row.year << "," << row.catch_mt << ","
              << row.biomass << "," << row.index << "\n";
    ++printed;
  }
}

inline void print_opakapaka_output_manifest() {
  std::cout << "\nWrote outputs:\n";
  std::cout << "  examples/NMFS/pifsc_opakapaka/outputs/"
               "synthetic_fit_summary.csv\n";
  std::cout << "  examples/NMFS/pifsc_opakapaka/outputs/"
               "synthetic_projection_scenarios.csv\n";
}

} // namespace opakapaka_example

using opakapaka_example::print_opakapaka_banner;
using opakapaka_example::print_opakapaka_fit_diagnostics;
using opakapaka_example::print_opakapaka_optimizer_structure;
using opakapaka_example::print_opakapaka_output_manifest;
using opakapaka_example::print_opakapaka_projection_preview;
