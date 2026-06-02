#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250}"
Rscript examples/tmb_state_space_surplus/benchmark_scaled_fixed_theta_tmb.R "$REPS" "$LENGTHS"
