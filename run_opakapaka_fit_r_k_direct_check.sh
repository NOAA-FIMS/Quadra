#!/usr/bin/env bash
set -euo pipefail

./inspect_opakapaka_fit_r_k_direct.sh

echo
echo "== O3 build after fitting log_r/log_K =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka_fit_r_k_check

echo
echo "== Run Opakapaka fit r/K check binary =="
./build/examples/pifsc_opakapaka_fit_r_k_check
