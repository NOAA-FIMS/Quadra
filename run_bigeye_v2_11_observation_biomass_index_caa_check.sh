#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/11_observation_biomass_index_caa/bigeye_v2_11_observation_biomass_index_caa_check.cpp \
  -o build/examples/bigeye_v2_11_observation_biomass_index_caa_check

./build/examples/bigeye_v2_11_observation_biomass_index_caa_check
