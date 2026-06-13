#include "red_snapper_model.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    out.push_back(item);
  }
  return out;
}

std::vector<sefsc_red_snapper::Observation> read_observations(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Could not open observations CSV: " + path);
  }

  std::string line;
  std::getline(in, line);  // header

  std::vector<sefsc_red_snapper::Observation> out;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const auto fields = split_csv_line(line);
    if (fields.size() != 13) {
      throw std::runtime_error("Expected 13 columns in observations CSV");
    }

    sefsc_red_snapper::Observation obs;
    obs.year = std::stoi(fields[0]);
    obs.catch_mt = std::stod(fields[1]);
    obs.index = std::stod(fields[2]);
    for (std::size_t a = 0; a < obs.age_comp.size(); ++a) {
      obs.age_comp[a] = std::stod(fields[3 + a]);
    }

    double age_comp_sum = 0.0;
    for (double v : obs.age_comp) {
      age_comp_sum += v;
    }
    if (age_comp_sum > 0.0) {
      for (double &v : obs.age_comp) {
        v /= age_comp_sum;
      }
    }
    out.push_back(obs);
  }
  return out;
}

void write_derived_quantities(
    const std::string& path,
    const std::vector<sefsc_red_snapper::DerivedRow>& rows) {
  std::ofstream out(path);
  out << "year,biomass,ssb_proxy,depletion,F_proxy,index_hat\n";
  out << std::fixed << std::setprecision(6);
  for (const auto& row : rows) {
    out << row.year << "," << row.biomass << "," << row.ssb_proxy << ","
        << row.depletion << "," << row.f_proxy << "," << row.index_hat << "\n";
  }
}

}  // namespace

int main() {
  const std::string input_path =
      "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv";
  const std::string output_path =
      "examples/NMFS/sefsc_red_snapper/outputs/level0_derived_quantities.csv";

  auto observations = read_observations(input_path);
  sefsc_red_snapper::RedSnapperModel model(observations);

  // Fixed placeholder values. Next patch should estimate these.
  const double log_r0 = std::log(1400.0);
  const double log_q = std::log(0.001);
  const double log_f = std::log(0.25);

  auto trajectory = model.deterministic_trajectory(log_r0, log_q, log_f);
  write_derived_quantities(output_path, trajectory);

  std::cout << "SEFSC red-snapper-style Level-0 scaffold\n";
  std::cout << "observations: " << observations.size() << "\n";
  std::cout << "wrote: " << output_path << "\n";

  if (!trajectory.empty()) {
    const auto& last = trajectory.back();
    std::cout << "terminal biomass: " << last.biomass << "\n";
    std::cout << "terminal depletion: " << last.depletion << "\n";
  }

  return 0;
}
