#!/usr/bin/env bash
set -euo pipefail

./inspect_opakapaka_empty_namespace_repair.sh

echo
echo "== O3 build after empty anonymous namespace repair =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka_refactor_check

echo
echo "== Run Opakapaka refactor check binary =="
./build/examples/pifsc_opakapaka_refactor_check
