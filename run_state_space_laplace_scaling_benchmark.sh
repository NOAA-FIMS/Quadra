#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
LENGTHS="${2:-1000,2000,5000,10000,20000}"

STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="benchmarks/state_space_laplace_scaling_${STAMP}"
mkdir -p "$OUTDIR"

RAW="${OUTDIR}/laplace_evaluator_scaling_raw.txt"
SUMMARY="${OUTDIR}/summary.tsv"

if [[ ! -x "./run_state_space_laplace_evaluator_benchmark.sh" ]]; then
  echo "ERROR: run_state_space_laplace_evaluator_benchmark.sh not found or not executable"
  exit 1
fi

cat > "${OUTDIR}/README.txt" <<README
State-space LaplaceEvaluator scaling benchmark
==============================================

Generated: ${STAMP}
Reps:      ${REPS}
Lengths:   ${LENGTHS}

Purpose:
  Test whether the optimized Quadra state-space Laplace runtime scales
  approximately linearly with the number of latent states.

Columns in summary.tsv:
  n
  objective
  lbfgs_ms
  cold_newton_ms
  warm_newton_ms
  warm_speedup_vs_lbfgs
  cold_iterations
  warm_iterations
  grad_norm
  objective_diff
README

echo "== State-space LaplaceEvaluator scaling benchmark =="
echo "reps = ${REPS}"
echo "lengths = ${LENGTHS}"
echo

./run_state_space_laplace_evaluator_benchmark.sh "$REPS" "$LENGTHS" | tee "$RAW"

printf "n\tobjective\tlbfgs_ms\tcold_newton_ms\twarm_newton_ms\twarm_speedup_vs_lbfgs\tcold_iterations\twarm_iterations\tgrad_norm\tobjective_diff\n" > "$SUMMARY"

awk '
  /^[[:space:]]*[0-9]+[[:space:]]/ {
    # n objective lbfgs cold warm speedup cold_it warm_it grad diff
    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10
  }
' "$RAW" >> "$SUMMARY"

echo
echo "Scaling benchmark complete."
echo "Output directory:"
echo "  ${OUTDIR}"
echo
echo "Summary:"
echo "  ${SUMMARY}"
echo
column -t -s $'\t' "$SUMMARY" || cat "$SUMMARY"
