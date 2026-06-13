#include "red_snapper_objective.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void write_objective_summary(
    const std::string& path,
    const sefsc_red_snapper::ObjectiveBreakdown& obj,
    const sefsc_red_snapper::AgeStructuredParams& params) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open objective summary CSV: " + path);
  }

  out << "field,value\n";
  out << std::setprecision(12);
  out << "objective_total," << obj.total << "\n";
  out << "index_nll," << obj.index_nll << "\n";
  out << "catch_nll," << obj.catch_nll << "\n";
  out << "n_index," << obj.n_index << "\n";
  out << "n_catch," << obj.n_catch << "\n";
  out << "log_r0," << params.log_r0 << "\n";
  out << "r0," << std::exp(params.log_r0) << "\n";
  out << "log_m," << params.log_m << "\n";
  out << "m," << std::exp(params.log_m) << "\n";
  out << "log_fbar," << params.log_fbar << "\n";
  out << "fbar," << std::exp(params.log_fbar) << "\n";
  out << "log_q," << params.log_q << "\n";
  out << "q," << std::exp(params.log_q) << "\n";
  out << "sel_a50," << params.sel_a50 << "\n";
  out << "sel_slope," << params.sel_slope << "\n";
}

}  // namespace

int main() {
  const std::string input_path =
      "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv";
  const std::string summary_path =
      "examples/NMFS/sefsc_red_snapper/outputs/objective_summary.csv";

  const auto observations = sefsc_red_snapper::read_observations(input_path);

  sefsc_red_snapper::AgeStructuredParams params;
  sefsc_red_snapper::ObjectiveOptions options;

  const auto breakdown =
      sefsc_red_snapper::evaluate_objective_breakdown(observations, params,
                                                      options);

  write_objective_summary(summary_path, breakdown, params);

  std::cout << "SEFSC red-snapper-style objective scaffold\n";
  std::cout << "objective_total: " << breakdown.total << "\n";
  std::cout << "index_nll:       " << breakdown.index_nll << "\n";
  std::cout << "catch_nll:       " << breakdown.catch_nll << "\n";
  std::cout << "wrote:           " << summary_path << "\n";

  return 0;
}
