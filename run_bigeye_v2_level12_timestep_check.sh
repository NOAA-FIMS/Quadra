#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level12_timestep_check/bigeye_v2_level12_timestep_check.cpp \
  -o build/examples/bigeye_v2_level12_timestep_check

./build/examples/bigeye_v2_level12_timestep_check
