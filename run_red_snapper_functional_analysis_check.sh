#!/usr/bin/env bash
set -euo pipefail

./inspect_red_snapper_functional_analysis_patch.sh

echo
echo "== O3 build after Red Snapper functional analysis =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp \
  -o build/examples/sefsc_red_snapper_functional_analysis_check

echo
echo "== Run Red Snapper functional analysis check binary =="
./build/examples/sefsc_red_snapper_functional_analysis_check

echo
echo "== Functional analysis outputs =="
ls -lh \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_functional_analysis_report.txt \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_functional_analysis_report.csv

echo
echo "== Functional analysis text preview =="
sed -n '1,120p' \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_functional_analysis_report.txt

echo
echo "== Functional analysis CSV preview =="
sed -n '1,40p' \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_functional_analysis_report.csv
