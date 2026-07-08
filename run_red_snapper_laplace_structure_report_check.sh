#!/usr/bin/env bash
set -euo pipefail

./inspect_red_snapper_laplace_structure_report.sh

echo
echo "== O3 build after Red Snapper Laplace structure report =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp \
  -o build/examples/sefsc_red_snapper_laplace_structure_check

echo
echo "== Run Red Snapper Laplace structure report check binary =="
./build/examples/sefsc_red_snapper_laplace_structure_check

echo
echo "== Laplace structure outputs =="
ls -lh \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_laplace_structure_report.txt \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_laplace_structure_report.csv

echo
echo "== Laplace structure text preview =="
sed -n '1,120p' \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_laplace_structure_report.txt

echo
echo "== Laplace structure CSV preview =="
sed -n '1,40p' \
  examples/NMFS/sefsc_red_snapper/outputs/red_snapper_laplace_structure_report.csv
