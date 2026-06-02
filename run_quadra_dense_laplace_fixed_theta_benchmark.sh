#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/benchmark_dense_laplace_fixed_theta.cpp \
  -o build/examples/benchmark_dense_laplace_fixed_theta

./build/examples/benchmark_dense_laplace_fixed_theta "$REPS"
