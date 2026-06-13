#include "red_snapper_age_structured.hpp"

#include <iostream>

int main() {
  const std::string input_path =
      "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv";
  const std::string output_path =
      "examples/NMFS/sefsc_red_snapper/outputs/age_structured_deterministic_trajectory.csv";

  const auto observations = sefsc_red_snapper::read_observations(input_path);

  sefsc_red_snapper::AgeStructuredParams params;
  const auto rows =
      sefsc_red_snapper::run_deterministic_age_structured_model(observations,
                                                                params);

  sefsc_red_snapper::write_age_structured_rows(output_path, rows);

  std::cout << "SEFSC red-snapper-style deterministic age-structured model\n";
  std::cout << "observations: " << observations.size() << "\n";
  std::cout << "wrote: " << output_path << "\n";

  if (!rows.empty()) {
    const auto& terminal = rows.back();
    std::cout << "terminal total biomass: " << terminal.total_biomass << "\n";
    std::cout << "terminal SSB proxy: " << terminal.ssb_proxy << "\n";
    std::cout << "terminal depletion: " << terminal.depletion << "\n";
  }

  return 0;
}
