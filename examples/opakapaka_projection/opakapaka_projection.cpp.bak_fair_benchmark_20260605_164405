#include "opakapaka_model.hpp"

#include <iomanip>
#include <iostream>
#include <chrono>

int main() {
  using namespace opakapaka_example;

  std::cout << "Synthetic opakapaka-style fit + projection example\n";
  std::cout << "==================================================\n\n";
  std::cout << "Synthetic and public-data-safe. Not an official assessment.\n\n";

  auto data = make_synthetic_opakapaka_data();

  OpakapakaProjectionModel model(data);
  auto params = model.initial_parameters();

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  // Public Quadra workflow:
  //   instantiate model -> optimize_lbfgs -> inspect fit -> project
  const auto fit_start = std::chrono::steady_clock::now();
  auto fit = quadra::optimize_lbfgs(model, params, opts);
  const auto fit_stop = std::chrono::steady_clock::now();
  const double fit_runtime_ms =
      std::chrono::duration<double, std::milli>(fit_stop - fit_start).count();

  ProjectionOptions projection_options;
  projection_options.start_year = data.back().year + 1;
  projection_options.years = 10;
  projection_options.scenarios = {
      {"zero_catch", 0.0},
      {"status_quo", 1.0},
      {"low_catch", 0.75},
      {"high_catch", 1.25},
  };

  auto projection = model.project(fit, projection_options);

  std::cout << "\nFit diagnostics\n";
  std::cout << "---------------\n";
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "objective          " << fit.value << "\n";
  std::cout << "grad_norm          " << fit.grad_norm << "\n";
  std::cout << "runtime_ms         " << fit_runtime_ms << "\n";
  std::cout << "iterations         " << fit.iterations << "\n";
  std::cout << "converged          " << (fit.converged ? "yes" : "no")
            << "\n";
  std::cout << "message            " << fit.message << "\n";
  std::cout << "log_q              " << fit.par.at(0) << "\n";
  std::cout << "q                  " << std::exp(fit.par.at(0)) << "\n";

  std::cout << "\nOptimizer structure diagnostics\n";
  std::cout << "-------------------------------\n";
  std::cout << "random effects     " << fit.pattern.random_effect_count << "\n";
  std::cout << "pattern available  " << (fit.pattern.available ? "yes" : "no")
            << "\n";
  std::cout << "detected structure " << fit.pattern.detected_structure << "\n";
  std::cout << "Laplace backend    " << fit.pattern.backend << "\n";
  std::cout << "random solver      " << fit.pattern.solver << "\n";
  std::cout << "complexity         " << fit.pattern.complexity << "\n";
  std::cout << "bandwidth          " << fit.pattern.bandwidth << "\n";
  std::cout << "Hessian nonzeros   " << fit.pattern.nonzeros << "\n";

  std::cout << "\nProjection preview\n";
  std::cout << "------------------\n";
  std::cout << "scenario,year,catch_mt,biomass,index\n";
  int printed = 0;
  for (const auto &row : projection) {
    if (printed >= 12) {
      break;
    }
    std::cout << row.scenario << "," << row.year << "," << row.catch_mt
              << "," << row.biomass << "," << row.index << "\n";
    ++printed;
  }

  write_fit_summary_csv(
      "examples/opakapaka_projection/outputs/synthetic_fit_summary.csv", fit);
  write_projection_csv(
      "examples/opakapaka_projection/outputs/synthetic_projection_scenarios.csv",
      projection);

  std::cout << "\nWrote outputs:\n";
  std::cout << "  examples/opakapaka_projection/outputs/"
               "synthetic_fit_summary.csv\n";
  std::cout << "  examples/opakapaka_projection/outputs/"
               "synthetic_projection_scenarios.csv\n";

  return 0;
}
