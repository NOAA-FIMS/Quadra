#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level01_population_check/bigeye_v2_level01_population_check.cpp \
  -o build/examples/bigeye_v2_level01_population_check

./build/examples/bigeye_v2_level01_population_check
