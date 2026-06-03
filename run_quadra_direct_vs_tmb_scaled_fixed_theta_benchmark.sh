#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"

mkdir -p build/examples

echo "== Quadra scaled direct latent-state tridiagonal =="
set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/benchmark_latent_tridiagonal_direct_scaled.cpp \
  -o build/examples/benchmark_latent_tridiagonal_direct_scaled
./build/examples/benchmark_latent_tridiagonal_direct_scaled "$REPS" "$LENGTHS"
set +x

echo
echo "== TMB scaled AD/Laplace =="
if [[ -x ./run_tmb_scaled_state_space_surplus_benchmark.sh ]]; then
  ./run_tmb_scaled_state_space_surplus_benchmark.sh "$REPS" "$LENGTHS"
else
  ./run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh "$REPS" "$LENGTHS" | sed -n '/== TMB scaled AD\/Laplace ==/,$p'
fi
