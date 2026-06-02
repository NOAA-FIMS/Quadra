#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"

echo "== Quadra latent-state tridiagonal Laplace =="
./run_quadra_latent_tridiagonal_laplace_benchmark.sh "$REPS"

echo
echo "== TMB AD/Laplace =="
./run_tmb_fixed_theta_benchmark.sh "$REPS"
