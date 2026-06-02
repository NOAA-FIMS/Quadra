#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250}"

echo "== Quadra scaled latent-state tridiagonal =="
./run_quadra_scaled_latent_tridiagonal_benchmark.sh "$REPS" "$LENGTHS"

echo
echo "== TMB scaled AD/Laplace =="
./run_tmb_scaled_fixed_theta_benchmark.sh "$REPS" "$LENGTHS"
