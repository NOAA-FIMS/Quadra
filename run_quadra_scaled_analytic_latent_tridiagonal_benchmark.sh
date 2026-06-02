#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250}"
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/benchmark_latent_tridiagonal_analytic_scaled.cpp \
  -o build/examples/benchmark_latent_tridiagonal_analytic_scaled

./build/examples/benchmark_latent_tridiagonal_analytic_scaled "$REPS" "$LENGTHS"
