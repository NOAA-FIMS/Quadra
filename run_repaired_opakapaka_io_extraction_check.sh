#!/usr/bin/env bash
set -euo pipefail

./inspect_repaired_opakapaka_io_extraction.sh

echo
echo "== O3 build after Opakapaka I/O repair =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_projection.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka_refactor_check

echo
echo "== Run Opakapaka refactor check binary =="
./build/examples/pifsc_opakapaka_refactor_check
