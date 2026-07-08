#!/usr/bin/env bash
set -euo pipefail

./inspect_opakapaka_reference_points.sh

echo
echo "== O3 build after Opakapaka reference points =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka_reference_points_check

echo
echo "== Run Opakapaka reference points check binary =="
./build/examples/pifsc_opakapaka_reference_points_check

echo
echo "== reference_points.csv =="
cat examples/NMFS/pifsc_opakapaka/outputs/reference_points.csv
