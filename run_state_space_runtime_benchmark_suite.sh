#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
LENGTHS="${2:-25,50,100,250,500,1000}"

STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="benchmarks/state_space_runtime_suite_${STAMP}"
mkdir -p "$OUTDIR"

SUMMARY="${OUTDIR}/summary.tsv"

cat > "${OUTDIR}/README.txt" <<README
State-space runtime benchmark suite
===================================

Generated: ${STAMP}
Reps:      ${REPS}
Lengths:   ${LENGTHS}

This suite compares the staged evolution of the Quadra state-space Laplace path.

Included when available:
1. analytic direct-runtime benchmark
2. direct runtime timing breakdown
3. tridiagonal Newton solver benchmark
4. warm-start Newton solver benchmark
5. runtime-owned random state benchmark
6. solver-policy PersistentLatentStateRuntime benchmark
7. LaplaceEvaluator benchmark
8. TMB benchmark

Interpretation:
- LBFGS timings are the baseline for random-effect optimization.
- cold Newton tests exact tridiagonal Newton from deterministic initial values.
- warm Newton tests cached xhat reuse.
- LaplaceEvaluator tests the reusable runtime API with policy-driven solve/evaluation.
README

printf "benchmark\tn\tobjective\tlbfgs_ms\tcold_newton_ms\twarm_newton_ms\tspeedup\tgrad_norm\tobjective_diff\n" > "$SUMMARY"

run_and_capture() {
  local label="$1"
  local script="$2"
  local outfile="$3"

  if [[ -x "$script" ]]; then
    echo
    echo "== ${label} =="
    "./${script}" "$REPS" "$LENGTHS" | tee "${OUTDIR}/${outfile}"
  else
    echo
    echo "== ${label} =="
    echo "SKIP: ${script} not found or not executable" | tee "${OUTDIR}/${outfile}"
  fi
}

append_warm_table_summary() {
  local label="$1"
  local file="$2"

  if [[ ! -f "$file" ]]; then
    return
  fi

  awk -v label="$label" '
    /^[[:space:]]*[0-9]+[[:space:]]/ {
      n=$1;
      objective=$2;
      lbfgs=$3;
      cold=$4;
      warm=$5;
      speedup=$6;
      grad=$9;
      diff=$10;
      print label "\t" n "\t" objective "\t" lbfgs "\t" cold "\t" warm "\t" speedup "\t" grad "\t" diff;
    }
  ' "$file" >> "$SUMMARY"
}

run_and_capture "Analytic direct-runtime benchmark" "run_state_space_analytic_direct_runtime_scaled_benchmark.sh" "analytic_direct_runtime.txt"

run_and_capture "Direct runtime timing breakdown" "run_state_space_direct_runtime_breakdown_benchmark.sh" "direct_runtime_breakdown.txt"

run_and_capture "Tridiagonal Newton solver benchmark" "run_state_space_tridiagonal_newton_solver_benchmark.sh" "tridiagonal_newton_solver.txt"

run_and_capture "Warm-start Newton solver benchmark" "run_state_space_warmstart_newton_solver_benchmark.sh" "warmstart_newton_solver.txt"
append_warm_table_summary "warmstart_newton" "${OUTDIR}/warmstart_newton_solver.txt"

run_and_capture "Runtime-owned random state benchmark" "run_state_space_runtime_owned_random_state_benchmark.sh" "runtime_owned_random_state.txt"
append_warm_table_summary "runtime_owned_random_state" "${OUTDIR}/runtime_owned_random_state.txt"

run_and_capture "Solver-policy runtime benchmark" "run_state_space_solver_policy_runtime_benchmark.sh" "solver_policy_runtime.txt"
append_warm_table_summary "solver_policy_runtime" "${OUTDIR}/solver_policy_runtime.txt"

run_and_capture "LaplaceEvaluator benchmark" "run_state_space_laplace_evaluator_benchmark.sh" "laplace_evaluator.txt"
append_warm_table_summary "laplace_evaluator" "${OUTDIR}/laplace_evaluator.txt"

if [[ -x "./run_tmb_scaled_fixed_theta_benchmark.sh" ]]; then
  echo
  echo "== TMB scaled fixed-theta benchmark =="
  ./run_tmb_scaled_fixed_theta_benchmark.sh "$REPS" "$LENGTHS" | tee "${OUTDIR}/tmb_scaled_fixed_theta.txt"
else
  echo
  echo "== TMB scaled fixed-theta benchmark =="
  echo "SKIP: run_tmb_scaled_fixed_theta_benchmark.sh not found or not executable" | tee "${OUTDIR}/tmb_scaled_fixed_theta.txt"
fi

echo
echo "Benchmark suite complete."
echo "Output directory:"
echo "  ${OUTDIR}"
echo
echo "Summary:"
echo "  ${SUMMARY}"
echo
echo "Quick view:"
column -t -s $'\t' "$SUMMARY" || cat "$SUMMARY"
