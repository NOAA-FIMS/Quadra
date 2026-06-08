#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-3}"
LENGTHS="${2:-25,50}"
AGES="${3:-10}"

echo "== Quadra dense Laplace validation =="
./run_quadra_age_structured_dense_laplace_validation.sh "$REPS" "$LENGTHS" "$AGES"

echo
echo "== TMB AD/Laplace =="
./run_tmb_age_structured_recruitment_benchmark.sh "$REPS" "$LENGTHS" "$AGES"
