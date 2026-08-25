#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples
"${CXX:-c++}" -std=c++17 -O2 -I. -Iexternal/eigen \
  examples/surplus_production/run_surplus_production.cpp \
  -o build/examples/run_surplus_production
./build/examples/run_surplus_production
