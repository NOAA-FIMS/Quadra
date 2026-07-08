#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level09_optimizer_check/bigeye_v2_level09_optimizer_check.cpp \
  -o build/examples/bigeye_v2_level09_optimizer_check

./build/examples/bigeye_v2_level09_optimizer_check
