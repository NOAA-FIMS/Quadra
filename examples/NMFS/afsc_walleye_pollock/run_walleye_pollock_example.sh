#!/usr/bin/env bash
set -euo pipefail
mkdir -p build/examples
clang++ -std=c++17 -O3 -I"external/eigen/" \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock_adgraph_global.cpp \
  -o build/examples/afsc_walleye_pollock
./build/examples/afsc_walleye_pollock
