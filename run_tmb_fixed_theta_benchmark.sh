#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
Rscript examples/tmb_state_space_surplus/benchmark_fixed_theta_tmb.R "$REPS"
