#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/16_movement_package_caa/bigeye_v2_16_movement_package_caa_check.cpp \
  -o build/examples/bigeye_v2_16_movement_package_caa_check

./build/examples/bigeye_v2_16_movement_package_caa_check
