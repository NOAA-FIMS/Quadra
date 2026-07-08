#!/usr/bin/env bash
set -euo pipefail

./inspect_red_snapper_reference_points.sh

echo
echo "== O3 build after Red Snapper reference points =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp \
  -o build/examples/sefsc_red_snapper_reference_points_check

echo
echo "== Run Red Snapper reference points check binary =="
./build/examples/sefsc_red_snapper_reference_points_check

echo
echo "== Red Snapper reference points =="
cat examples/NMFS/sefsc_red_snapper/outputs/red_snapper_reference_points.csv
