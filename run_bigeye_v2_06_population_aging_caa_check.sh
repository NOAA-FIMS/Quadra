#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/06_population_aging_caa/bigeye_v2_06_population_aging_caa_check.cpp \
  -o build/examples/bigeye_v2_06_population_aging_caa_check

./build/examples/bigeye_v2_06_population_aging_caa_check
