#include "../diagnostics/bigeye_functional_analysis_diagnostics.hpp"
#include "../objective/bigeye_quadra_objective.hpp"
#include "../reports/bigeye_report_suite.hpp"
#include "bigeye_age_structured.hpp"

#include "../../../../../core/optimizer.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
  const std::string input_path =
      "examples/NMFS/pifsc_bigeye_tuna/level0_single_region/data/"
      "synthetic_bigeye_level0_observations.csv";
  const auto report_paths = pifsc_bigeye_tuna::default_bigeye_report_paths();
  const auto observations = pifsc_bigeye_tuna::read_observations(input_path);

  pifsc_bigeye_tuna::BigeyeQuadraObjective objective(observations);

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

  pifsc_bigeye_tuna::write_bigeye_report_suite(report_paths, observations,
                                               objective, params, fit);
  pifsc_bigeye_tuna::write_bigeye_functional_analysis_report(
      "examples/NMFS/pifsc_bigeye_tuna/level0_single_region/outputs/"
      "bigeye_functional_analysis_report.txt",
      "examples/NMFS/pifsc_bigeye_tuna/level0_single_region/outputs/"
      "bigeye_functional_analysis_report.csv",
      objective, params, fit);
  std::cout
      << "PIFSC bigeye-tuna-style Quadra Laplace recruitment-deviation fit\n";
  std::cout << "objective:  " << fit.value << "\n";
  std::cout << "grad_norm:  " << fit.grad_norm << "\n";
  std::cout << "converged:  " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message:    " << fit.message << "\n";
  std::cout << "wrote:      " << report_paths.summary << "\n";
  std::cout << "wrote:      " << report_paths.trajectory << "\n";
  std::cout << "wrote:      " << report_paths.residual_diagnostics << "\n";
  std::cout << "wrote:      " << report_paths.selectivity << "\n";
  std::cout << "wrote:      " << report_paths.recruitment_deviations << "\n";
  std::cout << "wrote:      " << report_paths.objective_components << "\n";
  std::cout << "wrote:      " << report_paths.laplace_structure_text << "\n";
  std::cout << "wrote:      " << report_paths.laplace_structure_csv << "\n";
  std::cout << "wrote:      "
            << "examples/NMFS/pifsc_bigeye_tuna/level0_single_region/outputs/"
               "bigeye_functional_analysis_report.txt\n";

  std::cout << "wrote:      "
            << "examples/NMFS/pifsc_bigeye_tuna/level0_single_region/outputs/"
               "bigeye_functional_analysis_report.csv\n";
  std::cout << "wrote:      " << report_paths.reference_points << "\n";
  return 0;
}
