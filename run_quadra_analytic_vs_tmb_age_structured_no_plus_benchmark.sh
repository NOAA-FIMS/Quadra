#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

echo "== Quadra no-plus age-structured analytic banded Laplace =="
./run_quadra_age_structured_no_plus_analytic_banded_benchmark.sh "$REPS" "$LENGTHS" "$AGES"

echo
echo "== TMB no-plus age-structured AD/Laplace =="
./run_tmb_age_structured_no_plus_benchmark.sh "$REPS" "$LENGTHS" "$AGES"
