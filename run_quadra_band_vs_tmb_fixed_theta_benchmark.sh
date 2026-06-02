#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
BANDWIDTH="${2:-3}"

echo "== Quadra band-Huu finite-difference Laplace =="
./run_quadra_band_laplace_fixed_theta_benchmark.sh "$REPS" "$BANDWIDTH"

echo
echo "== TMB AD/Laplace =="
./run_tmb_fixed_theta_benchmark.sh "$REPS"
