#!/usr/bin/env bash
set -euo pipefail

./inspect_red_snapper_fit_reports_extraction.sh

echo
echo "== O3 build after Red Snapper report extraction =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp \
  -o build/examples/sefsc_red_snapper_fit_reports_check

echo
echo "== Run Red Snapper fit report extraction check binary =="
./build/examples/sefsc_red_snapper_fit_reports_check
