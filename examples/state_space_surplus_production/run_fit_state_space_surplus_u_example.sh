#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples
"${CXX:-c++}" -std=c++17 -O2 -I. -Iexternal/eigen \
  -Iexternal/LBFGSpp/include \
  examples/state_space_surplus_production/fit_state_space_surplus_u.cpp \
  -o build/examples/fit_state_space_surplus_u
./build/examples/fit_state_space_surplus_u
