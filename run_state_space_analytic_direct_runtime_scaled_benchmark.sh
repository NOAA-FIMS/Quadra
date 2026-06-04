#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"

mkdir -p build/examples benchmarks

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -I. \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/benchmark_latent_tridiagonal_analytic_direct_runtime_scaled.cpp \
  -o build/examples/benchmark_latent_tridiagonal_analytic_direct_runtime_scaled

./build/examples/benchmark_latent_tridiagonal_analytic_direct_runtime_scaled "$REPS" "$LENGTHS" \
  | tee "benchmarks/state_space_analytic_direct_runtime_scaled_reps${REPS}_n${LENGTHS//,/}.txt"
