#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"

echo "== Quadra dense finite-difference Laplace =="
./run_quadra_dense_laplace_fixed_theta_benchmark.sh "$REPS"

echo
echo "== TMB AD/Laplace =="
./run_tmb_fixed_theta_benchmark.sh "$REPS"
