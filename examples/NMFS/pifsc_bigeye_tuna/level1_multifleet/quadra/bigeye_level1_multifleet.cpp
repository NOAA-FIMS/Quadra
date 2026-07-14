#include "../diagnostics/bigeye_fixed_effect_geometry.hpp"
#include "../diagnostics/bigeye_functional_analysis_diagnostics.hpp"
#include "../objective/bigeye_quadra_objective.hpp"
#include "../reports/bigeye_report_suite.hpp"
#include "bigeye_age_structured.hpp"

#include "../../../../../core/optimizer.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<pifsc_bigeye_tuna::Observation>
read_multifleet_observations_as_aggregate(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Could not open multifleet observations CSV: " +
                             path);
  }

  std::string line;
  std::getline(in, line);

  struct Accumulator {
    double catch_mt = 0.0;
    double index_sum = 0.0;
    int index_count = 0;
  };

  std::map<int, Accumulator> by_year;

  while (std::getline(in, line)) {
    if (line.empty())
      continue;

    std::stringstream ss(line);
    std::string year_s;
    std::string fleet;
    std::string catch_s;
    std::string index_s;

    std::getline(ss, year_s, ',');
    std::getline(ss, fleet, ',');
    std::getline(ss, catch_s, ',');
    std::getline(ss, index_s, ',');

    const int year = std::stoi(year_s);
    auto &acc = by_year[year];
    acc.catch_mt += std::stod(catch_s);
    acc.index_sum += std::stod(index_s);
    acc.index_count += 1;
  }

  std::vector<pifsc_bigeye_tuna::Observation> out;
  out.reserve(by_year.size());

  for (const auto &kv : by_year) {
    const int year = kv.first;
    const auto &acc = kv.second;
    const double index =
        (acc.index_count > 0)
            ? acc.index_sum / static_cast<double>(acc.index_count)
            : 0.0;

    std::array<double, pifsc_bigeye_tuna::kAges> age_comp{};
    const double inv_n_age =
        1.0 / static_cast<double>(pifsc_bigeye_tuna::kAges);
    for (double &p_age : age_comp) {
      p_age = inv_n_age;
    }

    out.push_back({year, acc.catch_mt, index, age_comp});
  }

  return out;
}

} // namespace

int main() {
  const std::string input_path =
      "examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/data/"
      "synthetic_bigeye_level1_multifleet_observations.csv";
  const auto report_paths = pifsc_bigeye_tuna::default_bigeye_report_paths();
  const auto observations =
      read_multifleet_observations_as_aggregate(input_path);

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
      "examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/"
      "bigeye_level1_functional_analysis_report.txt",
      "examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/"
      "bigeye_level1_functional_analysis_report.csv",
      objective, params, fit);
  pifsc_bigeye_tuna::write_bigeye_fixed_effect_geometry_report(
      "examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/"
      "bigeye_level1_fixed_effect_geometry_report.txt",
      "examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/"
      "bigeye_level1_fixed_effect_geometry_report.csv",
      objective, params, fit, opts);
  std::cout << "PIFSC bigeye-tuna-style Level 1 multi-fleet Quadra Laplace "
               "recruitment-deviation fit\n";
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
            << "examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/"
               "bigeye_level1_functional_analysis_report.txt\n";

  std::cout << "wrote:      "
            << "examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/"
               "bigeye_level1_functional_analysis_report.csv\n";
  std::cout << "wrote:      " << report_paths.reference_points << "\n";
  return 0;
}
