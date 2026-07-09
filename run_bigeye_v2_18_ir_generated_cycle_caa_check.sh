#!/usr/bin/env bash
set -euo pipefail

./generate_bigeye_v2_caa_assessment_cycle_from_ir.sh

mkdir -p build/examples

cat > build/examples/bigeye_v2_18_ir_generated_cycle_caa_check.cpp <<'CPP'
#include "examples/NMFS/pifsc_bigeye_tuna/v2/architecture/assessment/generated_assessment_cycle.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {
bool nearly_equal(double a, double b, double tol = 1.0e-6) {
  return std::abs(a - b) <= tol;
}
}

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 3;

  AssessmentParameters<double> parameters;
  parameters.life.log_m_young_offset = std::log(0.75);
  parameters.life.log_m_old_offset = std::log(0.65);

  parameters.populations.resize(1);
  parameters.populations[0].r0 = 1000.0;

  parameters.fleets.resize(1);
  parameters.fleets[0].sel_a50 = 5.0;
  parameters.fleets[0].sel_slope = 1.0;
  parameters.fleets[0].fbar = 0.2;
  parameters.fleets[0].q_index = 0.01;
  parameters.fleets[0].catch_sigma = 0.1;
  parameters.fleets[0].index_sigma = 0.2;

  AssessmentState<double> state;
  state.populations.resize(1);
  state.fleets.resize(1);

  state.populations[0].numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    state.populations[0].numbers_at_age[0][a] = 1000.0;
  }

  GeneratedAssessmentCycleFromIR{}(data, parameters, state);

  if (!nearly_equal(state.life.m_at_age[0], 0.3375)) {
    std::cerr << "FAIL: generated IR cycle life history\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.populations[0].numbers_at_age[1][1],
                    710.98976678322833)) {
    std::cerr << "FAIL: generated IR cycle population\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(state.fleets[0].z_at_age[4], 0.55)) {
    std::cerr << "FAIL: generated IR cycle fleet mortality\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA IR-generated AssessmentCycle regression\n";
  return EXIT_SUCCESS;
}
CPP

c++ -std=c++17 -O3 \
  -I. \
  build/examples/bigeye_v2_18_ir_generated_cycle_caa_check.cpp \
  -o build/examples/bigeye_v2_18_ir_generated_cycle_caa_check

./build/examples/bigeye_v2_18_ir_generated_cycle_caa_check
