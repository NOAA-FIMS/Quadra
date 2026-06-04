#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  -Iexamples/state_space_surplus_production \
  examples/state_space_surplus_production/analyze_state_space_hessian.cpp \
  -o build/examples/analyze_state_space_hessian

./build/examples/analyze_state_space_hessian
