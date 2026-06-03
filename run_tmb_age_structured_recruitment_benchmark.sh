#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

Rscript examples/tmb_age_structured_recruitment/benchmark_age_structured_recruitment_tmb.R "$REPS" "$LENGTHS" "$AGES"
