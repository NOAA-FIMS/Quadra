#include "red_snapper_age_structured.hpp"
#include "../objective/red_snapper_quadra_objective.hpp"
#include "../reports/red_snapper_fit_reports.hpp"

#include "../../../../core/optimizer.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
  const std::string input_path = "examples/NMFS/sefsc_red_snapper/data/"
                                 "synthetic_red_snapper_observations.csv";
  const std::string summary_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv";
  const std::string trajectory_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fitted_trajectory.csv";
  const std::string residual_diagnostics_path =
      "examples/NMFS/sefsc_red_snapper/outputs/"
      "quadra_fit_residual_diagnostics.csv";
  const std::string selectivity_path =
      "examples/NMFS/sefsc_red_snapper/outputs/selectivity_at_age.csv";
  const std::string recruitment_deviations_path =
      "examples/NMFS/sefsc_red_snapper/outputs/recruitment_deviations.csv";
  const std::string objective_components_path =
      "examples/NMFS/sefsc_red_snapper/outputs/"
      "quadra_fit_objective_components.csv";
  const auto observations = sefsc_red_snapper::read_observations(input_path);

  sefsc_red_snapper::RedSnapperQuadraObjective objective(observations);

  quadra::ParameterVector params;
  params.add({"log_r0", std::log(1200.0), quadra::ParameterTransform::Identity,
              false});
  params.add({"log_fbar", std::log(0.025), quadra::ParameterTransform::Identity,
              false});
  params.add({"log_q", std::log(0.00005), quadra::ParameterTransform::Identity,
              false});
  params.add(
      {"logit_sel_a50", 0.0, quadra::ParameterTransform::Identity, false});
  params.add({"log_sel_slope", std::log(1.2),
              quadra::ParameterTransform::Identity, false});

  for (std::size_t t = 0; t < observations.size(); ++t) {
    params.add({"log_rec_dev_" + std::to_string(t + 1), 0.0,
                quadra::ParameterTransform::Identity, true});
  }

  quadra::LaplaceOptions opts;

  auto fit = quadra::optimize_lbfgs(objective, params, opts);

  sefsc_red_snapper::write_fit_summary(summary_path, fit);
  write_fitted_trajectory(trajectory_path, observations, fit);
  write_residual_diagnostics(residual_diagnostics_path, observations, fit);
  write_selectivity_at_age(selectivity_path, fit);
  write_recruitment_deviations(recruitment_deviations_path, fit);
  write_objective_components(objective_components_path, observations, fit);
  std::cout
      << "SEFSC red-snapper-style Quadra Laplace recruitment-deviation fit\n";
  std::cout << "objective:  " << fit.value << "\n";
  std::cout << "grad_norm:  " << fit.grad_norm << "\n";
  std::cout << "converged:  " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message:    " << fit.message << "\n";
  std::cout << "wrote:      " << summary_path << "\n";
  std::cout << "wrote:      " << trajectory_path << "\n";
  std::cout << "wrote:      " << residual_diagnostics_path << "\n";
  std::cout << "wrote:      " << selectivity_path << "\n";
  std::cout << "wrote:      " << recruitment_deviations_path << "\n";
  std::cout << "wrote:      " << objective_components_path << "\n";
  return 0;
}
