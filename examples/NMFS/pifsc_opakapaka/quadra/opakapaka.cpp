#include "../../../../core/uncertainty/reporting.hpp"
#include "../../../../core/uncertainty/selected_inverse_diagonal.hpp"
#include "../data/opakapaka_io.hpp"
#include "../diagnostics/opakapaka_biomass_covariance_diagnostics.hpp"
#include "../diagnostics/opakapaka_logq_diagnostics.hpp"
#include "../diagnostics/opakapaka_projection_uncertainty.hpp"
#include "../diagnostics/opakapaka_random_effect_diagnostics.hpp"
#include "../optimization/opakapaka_logq_optimization.hpp"
#include "../reports/opakapaka_report_suite.hpp"
#include "drivers/opakapaka_driver_output.hpp"
#include "opakapaka_model.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
  using namespace opakapaka_example;

  std::cout << "Synthetic opakapaka-style fit + projection example\n";
  std::cout << "==================================================\n\n";
  std::cout
      << "Synthetic and public-data-safe. Not an official assessment.\n\n";

  auto data =
      read_opakapaka_history_csv("examples/NMFS/pifsc_opakapaka/data/"
                                 "synthetic_opakapaka_projection_data.csv");

  std::cout << "Loaded shared CSV fit rows: " << data.size() << "\n\n";

  OpakapakaProjectionModel model(data);
  auto params = model.initial_parameters();

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  // Public Quadra workflow:
  //   instantiate model -> optimize_lbfgs -> inspect fit -> project
  const auto fit_start = std::chrono::steady_clock::now();
  quadra::OptResult fit;
  bool primary_optimizer_converged = false;
  bool fallback_used = false;
  std::string primary_optimizer_name = "profiled scalar Laplace";
  std::string primary_optimizer_status = "not run";
  double primary_optimizer_grad_norm = std::numeric_limits<double>::quiet_NaN();

#ifndef OPAKAPAKA_USE_LBFGS_PRIMARY
  // Opakapaka has one fixed effect and twenty random effects. For this
  // geometry, the safeguarded profiled scalar Laplace optimizer is the
  // appropriate primary optimizer: it directly optimizes log_q while profiling
  // over the random effects and avoids quasi-Newton line-search pathologies.
  fit = quadra::optimize_lbfgs(model, params, opts);

  if (fit.converged) {
    fit.message = "converged with L-BFGS optimizer";
  }

  primary_optimizer_converged = fit.converged;
  primary_optimizer_status = fit.message;
  primary_optimizer_grad_norm = fit.grad_norm;
#else
  primary_optimizer_name = "L-BFGS";
  try {
    fit = quadra::optimize_lbfgs(model, params, opts);
    primary_optimizer_converged = fit.converged;
    primary_optimizer_status = fit.message;
    primary_optimizer_grad_norm = fit.grad_norm;
  } catch (const std::runtime_error &e) {
    const std::string msg = e.what();
    if (msg.find("line search") == std::string::npos &&
        msg.find("sufficiently decrease") == std::string::npos) {
      throw;
    }

    fallback_used = true;
    primary_optimizer_converged = false;
    primary_optimizer_status = msg;

    std::cout << "L-BFGS line-search stall detected in Opakapaka example. "
              << "Using local safeguarded one-dimensional log_q fallback.\n";

    fit = fit_log_q_fd_newton_fallback(model, params, opts,
                                       params.params.at(0).value);
  }
#endif

  const double fit_value_before_polish = fit.value;
  const double fit_grad_before_polish = fit.grad_norm;
  polish_single_logq_if_helpful(model, params, opts, fit);

  const bool polish_changed =
      std::abs(fit.value - fit_value_before_polish) > 1.0e-10 ||
      std::abs(fit.grad_norm - fit_grad_before_polish) > 1.0e-10;

#ifdef OPAKAPAKA_USE_LBFGS_PRIMARY
  fallback_used = fallback_used || polish_changed;
#else
  // In the default build, scalar optimization is primary. Optional scalar
  // polishing is still part of that primary scalar workflow, not a fallback.
  fallback_used = false;
  primary_optimizer_converged = fit.converged;
  primary_optimizer_status = fit.message;
  primary_optimizer_grad_norm = fit.grad_norm;
#endif

  const std::string convergence_status =
      primary_optimizer_converged && !fallback_used
          ? "primary_optimizer_converged"
          : (fallback_used ? "fallback_polished" : "not_converged");

  {
    std::ofstream state_out(
        "examples/NMFS/pifsc_opakapaka/outputs/quadra_fitted_states.csv");

    state_out << "index,log_B,B\n";

    for (std::size_t i = 0; i < fit.u_hat.size(); ++i) {
      state_out << i << "," << std::setprecision(15) << fit.u_hat[i] << ","
                << std::setprecision(15) << std::exp(fit.u_hat[i]) << "\n";
    }
  }

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

  const Eigen::SparseMatrix<double> Huu_final =
      compute_final_random_effect_hessian(model, params, opts, fit);
  const int final_hessian_nonzeros = static_cast<int>(Huu_final.nonZeros());

  print_opakapaka_fit_diagnostics(
      fit, fit_runtime_ms, convergence_status, primary_optimizer_name,
      fallback_used, primary_optimizer_converged, primary_optimizer_grad_norm,
      primary_optimizer_status);

  print_opakapaka_optimizer_structure(fit, final_hessian_nonzeros);
  print_opakapaka_projection_preview(projection);

  const auto final_h_uu =
      compute_final_random_effect_hessian(model, params, opts, fit);

  write_opakapaka_report_suite(model, params, opts, fit, data, projection,
                               final_h_uu);

  print_opakapaka_output_manifest();

  return 0;
}
