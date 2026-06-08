#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-3}"
LENGTHS="${2:-25,50}"
AGES="${3:-10}"
BANDWIDTH="${4:-9}"

echo "== Quadra banded vs dense Laplace validation =="
./run_quadra_age_structured_banded_laplace_validation.sh "$REPS" "$LENGTHS" "$AGES" "$BANDWIDTH"

echo
echo "== TMB AD/Laplace =="
./run_tmb_age_structured_recruitment_benchmark.sh "$REPS" "$LENGTHS" "$AGES"
