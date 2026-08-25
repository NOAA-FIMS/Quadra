#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples
"${CXX:-c++}" -std=c++17 -O3 -Iexternal/eigen \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka
./build/examples/pifsc_opakapaka
