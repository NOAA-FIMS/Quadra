#!/usr/bin/env bash
set -euo pipefail

./inspect_red_snapper_objective_extraction.sh

echo
echo "== O3 build after Red Snapper objective extraction =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp \
  -o build/examples/sefsc_red_snapper_objective_check

echo
echo "== Run Red Snapper objective extraction check binary =="
./build/examples/sefsc_red_snapper_objective_check
