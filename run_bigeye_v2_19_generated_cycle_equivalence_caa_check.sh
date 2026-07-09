#!/usr/bin/env bash
set -euo pipefail

./generate_bigeye_v2_caa_assessment_cycle_from_ir.sh

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/19_generated_cycle_equivalence_caa/bigeye_v2_19_generated_cycle_equivalence_caa_check.cpp \
  -o build/examples/bigeye_v2_19_generated_cycle_equivalence_caa_check

./build/examples/bigeye_v2_19_generated_cycle_equivalence_caa_check
