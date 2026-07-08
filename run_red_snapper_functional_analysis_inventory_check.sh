#!/usr/bin/env bash
set -euo pipefail

echo "== Current Red Snapper O3 build =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp \
  -o build/examples/sefsc_red_snapper_functional_inventory_check

echo
echo "== Current Red Snapper run =="
./build/examples/sefsc_red_snapper_functional_inventory_check
