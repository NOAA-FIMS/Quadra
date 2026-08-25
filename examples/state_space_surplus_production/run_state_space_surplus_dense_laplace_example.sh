#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples
"${CXX:-c++}" -std=c++17 -O2 -I. -Iexternal/eigen \
  -Iexternal/LBFGSpp/include \
  examples/state_space_surplus_production/laplace_state_space_surplus_dense.cpp \
  -o build/examples/laplace_state_space_surplus_dense
./build/examples/laplace_state_space_surplus_dense
