#include "../../core/laplace/model_analysis_report.hpp"

#include <Eigen/Sparse>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct ProjectionScenario {
  std::string name;
  int history_years = 0;
  int projection_years = 0;
  double quadra_warm_ms = 0.0;
  double quadra_rss_mb = 0.0;
};

Eigen::SparseMatrix<double>
make_first_order_markov_hessian(const int latent_states,
                                const double process_precision,
                                const double observation_precision) {
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(3 * latent_states));

  for (int i = 0; i < latent_states; ++i) {
    double diag = observation_precision + process_precision;

    if (i + 1 < latent_states) {
      diag += process_precision;
      const double off = -process_precision;
      triplets.emplace_back(i, i + 1, off);
      triplets.emplace_back(i + 1, i, off);
    }

    triplets.emplace_back(i, i, diag);
  }

  Eigen::SparseMatrix<double> H(latent_states, latent_states);
  H.setFromTriplets(triplets.begin(), triplets.end());
  H.makeCompressed();
  return H;
}

void print_projection_summary_row(
    const ProjectionScenario &scenario,
    const quadra::laplace::ModelAnalysisReport &report) {
  std::cout << scenario.name << "\t" << scenario.history_years << "\t"
            << scenario.projection_years << "\t" << report.random_effects
            << "\t" << quadra::laplace::ToString(report.structure) << "\t"
            << report.bandwidth << "\t"
            << quadra::laplace::ToString(report.solver) << "\t"
            << quadra::laplace::ToString(report.backend) << "\t"
            << report.complexity << "\t" << std::fixed << std::setprecision(3)
            << scenario.quadra_warm_ms << "\t" << scenario.quadra_rss_mb
            << "\n";
}

} // namespace

int main() {
  std::cout << "Opakapaka / Deep 7 projection structure demo\n";
  std::cout << "============================================\n\n";

  std::cout << "Synthetic, public-data-safe computational example.\n";
  std::cout << "Not an official stock assessment.\n\n";

  quadra::laplace::StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = false;
  opts.dense_size_cutoff = 0;
  opts.banded_width_cutoff = 64;

  const std::vector<ProjectionScenario> scenarios = {
      {"short_projection", 30, 30, 0.005613, 1.297},
      {"moderate_projection", 100, 100, 0.025100, 1.297},
      {"long_history_short_projection", 300, 100, 0.033646, 1.359},
      {"long_history_long_projection", 300, 300, 0.038004, 1.531},
      {"stress_2k", 1000, 1000, 0.135108, 1.938},
      {"stress_20k", 10000, 10000, 1.290508, 6.562},
      {"stress_100k", 50000, 50000, 6.645117, 57.609},
  };

  constexpr double process_precision = 1.0 / (0.15 * 0.15);
  constexpr double observation_precision = 1.0 / (0.20 * 0.20);

  std::cout << "Projection structure summary\n";
  std::cout << "scenario\thistory\tprojection\tlatent_"
               "states\tstructure\tbandwidth\tsolver\tbackend\tcomplexity\twarm"
               "_ms\trss_mb\n";

  for (const auto &scenario : scenarios) {
    const int latent_states =
        scenario.history_years + scenario.projection_years;
    const Eigen::SparseMatrix<double> H = make_first_order_markov_hessian(
        latent_states, process_precision, observation_precision);

    const auto report = quadra::laplace::analyze_hessian_structure(H, opts);
    print_projection_summary_row(scenario, report);
  }

  std::cout << "\nDetailed report for moderate_projection\n";
  std::cout << "---------------------------------------\n";

  const ProjectionScenario selected{"moderate_projection", 100, 100, 0.025100,
                                    1.297};
  const Eigen::SparseMatrix<double> H = make_first_order_markov_hessian(
      selected.history_years + selected.projection_years, process_precision,
      observation_precision);

  const auto report = quadra::laplace::analyze_hessian_structure(H, opts);
  report.Print(std::cout);
  std::cout << "\n";

  return 0;
}
