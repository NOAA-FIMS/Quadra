#!/usr/bin/env bash
set -euo pipefail

echo "== Add SEFSC red snapper synthetic data and initial model scaffold =="

BASE="examples/NMFS/sefsc_red_snapper"
mkdir -p "$BASE"/{data,quadra,tmb,outputs,validation}

cat > "$BASE/data/synthetic_red_snapper_observations.csv" <<'CSV'
year,catch_mt,index,age1,age2,age3,age4,age5,age6,age7,age8,age9,age10
1,220,0.82,0.18,0.21,0.19,0.15,0.10,0.07,0.04,0.03,0.02,0.01
2,230,0.86,0.17,0.22,0.19,0.15,0.10,0.07,0.04,0.03,0.02,0.01
3,245,0.89,0.16,0.22,0.20,0.15,0.10,0.07,0.04,0.03,0.02,0.01
4,260,0.91,0.15,0.21,0.21,0.16,0.10,0.07,0.04,0.03,0.02,0.01
5,275,0.93,0.15,0.20,0.21,0.16,0.11,0.07,0.04,0.03,0.02,0.01
6,290,0.95,0.14,0.20,0.21,0.17,0.11,0.07,0.04,0.03,0.02,0.01
7,305,0.96,0.14,0.19,0.21,0.17,0.11,0.08,0.04,0.03,0.02,0.01
8,315,0.94,0.13,0.19,0.21,0.17,0.12,0.08,0.04,0.03,0.02,0.01
9,320,0.91,0.13,0.18,0.21,0.18,0.12,0.08,0.04,0.03,0.02,0.01
10,330,0.88,0.12,0.18,0.21,0.18,0.12,0.08,0.05,0.03,0.02,0.01
11,335,0.84,0.12,0.17,0.21,0.18,0.13,0.08,0.05,0.03,0.02,0.01
12,340,0.81,0.11,0.17,0.20,0.19,0.13,0.09,0.05,0.03,0.02,0.01
13,330,0.80,0.12,0.17,0.20,0.18,0.13,0.09,0.05,0.03,0.02,0.01
14,320,0.82,0.13,0.18,0.20,0.18,0.12,0.09,0.05,0.03,0.02,0.01
15,310,0.85,0.14,0.18,0.20,0.17,0.12,0.08,0.05,0.03,0.02,0.01
16,300,0.89,0.15,0.19,0.20,0.17,0.11,0.08,0.05,0.03,0.02,0.01
17,295,0.93,0.16,0.19,0.20,0.16,0.11,0.08,0.05,0.03,0.02,0.01
18,285,0.97,0.17,0.20,0.19,0.16,0.10,0.08,0.05,0.03,0.02,0.01
19,275,1.01,0.18,0.20,0.19,0.15,0.10,0.08,0.05,0.03,0.01,0.01
20,265,1.04,0.19,0.21,0.18,0.15,0.10,0.07,0.05,0.03,0.01,0.01
CSV

cat > "$BASE/data/red_snapper_projection_scenarios.csv" <<'CSV'
scenario,projection_year,catch_mt
zero_catch,21,0
zero_catch,22,0
zero_catch,23,0
zero_catch,24,0
zero_catch,25,0
status_quo,21,265
status_quo,22,265
status_quo,23,265
status_quo,24,265
status_quo,25,265
high_catch,21,340
high_catch,22,340
high_catch,23,340
high_catch,24,340
high_catch,25,340
CSV

cat > "$BASE/quadra/red_snapper_model.hpp" <<'CPP'
#pragma once

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace sefsc_red_snapper {

struct Observation {
  int year = 0;
  double catch_mt = 0.0;
  double index = 0.0;
  std::array<double, 10> age_comp{};
};

struct ProjectionScenario {
  std::string scenario;
  int projection_year = 0;
  double catch_mt = 0.0;
};

struct DerivedRow {
  int year = 0;
  double biomass = 0.0;
  double ssb_proxy = 0.0;
  double depletion = 0.0;
  double f_proxy = 0.0;
  double index_hat = 0.0;
};

// Level-0 placeholder model:
// This is intentionally minimal. The next patch should replace this with
// Quadra AD/Laplace evaluation and recruitment deviations as random effects.
class RedSnapperModel {
 public:
  explicit RedSnapperModel(std::vector<Observation> obs)
      : observations_(std::move(obs)) {}

  const std::vector<Observation>& observations() const { return observations_; }

  std::vector<DerivedRow> deterministic_trajectory(double log_r0,
                                                    double log_q,
                                                    double log_f) const {
    const double r0 = std::exp(log_r0);
    const double q = std::exp(log_q);
    const double f = std::exp(log_f);

    std::vector<DerivedRow> out;
    out.reserve(observations_.size());

    double biomass = r0;
    const double unfished = r0;

    for (const auto& obs : observations_) {
      biomass = std::max(1.0, biomass + 0.25 * r0 - obs.catch_mt - 0.05 * biomass);
      DerivedRow row;
      row.year = obs.year;
      row.biomass = biomass;
      row.ssb_proxy = 0.35 * biomass;
      row.depletion = biomass / unfished;
      row.f_proxy = f * obs.catch_mt / std::max(1.0, biomass);
      row.index_hat = q * biomass;
      out.push_back(row);
    }

    return out;
  }

 private:
  std::vector<Observation> observations_;
};

}  // namespace sefsc_red_snapper
CPP

cat > "$BASE/quadra/red_snapper_level0.cpp" <<'CPP'
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
CPP

cat > "$BASE/tmb/red_snapper_tmb.cpp" <<'CPP'
// Placeholder TMB reference implementation for the SEFSC red-snapper-style example.
//
// The next milestone should implement the same likelihood and derived quantities
// as the Quadra model so objective values, estimates, random effects, and
// uncertainty outputs can be compared side by side.

template<class Type>
Type objective_function<Type>::operator()() {
  return Type(0.0);
}
CPP

cat > "$BASE/tmb/README.md" <<'MD'
# TMB Reference Implementation

This directory will contain the TMB comparison model for the SEFSC red-snapper-style example.

The current file is a placeholder and should not be used as a scientific reference yet.
MD

cat > "$BASE/validation/level0_checklist.md" <<'MD'
# Level-0 Checklist

- [ ] deterministic age-structured dynamics implemented
- [ ] synthetic catch observations read from `data/`
- [ ] synthetic index observations read from `data/`
- [ ] derived quantities written to `outputs/`
- [ ] minimal runner compiles from a clean checkout
- [ ] TMB reference implementation added
- [ ] Quadra/TMB comparison table added
MD

cat > "$BASE/outputs/.gitignore" <<'EOF'
*
!.gitignore
EOF

cat > "$BASE/run_red_snapper_level0.sh" <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/NMFS/sefsc_red_snapper/outputs

c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -o examples/NMFS/sefsc_red_snapper/quadra/red_snapper_level0 \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_level0.cpp

./examples/NMFS/sefsc_red_snapper/quadra/red_snapper_level0
SH
chmod +x "$BASE/run_red_snapper_level0.sh"

echo
echo "Created initial SEFSC red snapper scaffold."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_red_snapper_level0.sh"
echo
echo "Then inspect:"
echo "  head examples/NMFS/sefsc_red_snapper/outputs/level0_derived_quantities.csv"
