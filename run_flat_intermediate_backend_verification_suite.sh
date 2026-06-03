#!/usr/bin/env bash
set -euo pipefail

# run_flat_intermediate_backend_verification_suite.sh
#
# Verification suite for the experimental flat intermediate directional backend.
#
# Goal:
#   Confirm the large speedup did not silently change scalar-vs-batch Hdot,
#   workspace gradients, or exact-gradient behavior.

timestamp="$(date +%Y%m%d_%H%M%S)"
log="flat_intermediate_backend_verification_${timestamp}.log"

run_one() {
  local cmd="$1"
  echo
  echo "================================================================================"
  echo "RUN: $cmd"
  echo "================================================================================"
  bash -lc "$cmd"
}

{
  echo "Flat intermediate backend verification"
  echo "timestamp: $timestamp"
  echo

  run_one "./run_had_quadra_nonzero_batch_directional_test.sh"
  run_one "./run_had_quadra_directional_batch_propagation_test.sh"
  run_one "./run_had_graph_workspace_test.sh"
  run_one "./run_exact_gradient_workspace_test.sh"
  run_one "./run_laplace_exact_gradient_evaluator_test.sh"
  run_one "./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10"
  run_one "./run_sparse_rw1_flat_accumulator_trace_check.sh 10"

  echo
  echo "================================================================================"
  echo "Verification suite completed."
  echo "Log: $log"
  echo "================================================================================"
} 2>&1 | tee "$log"

echo
echo "Quick checks:"
grep -E "diff|grad diff|obj diff|trace diff|tests passed|diagnostic passed|Benchmark complete|flat intermediate:" "$log" || true

echo
echo "Saved log:"
echo "  $log"
