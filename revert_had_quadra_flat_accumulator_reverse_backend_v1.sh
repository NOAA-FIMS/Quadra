#!/usr/bin/env bash
set -euo pipefail

# revert_had_quadra_flat_accumulator_reverse_backend_v1.sh
#
# Reverts only the experimental flat-accumulator reverse backend changes in:
#   core/had_quadra.hpp
#
# It restores the newest backup created by:
#   install_had_quadra_flat_accumulator_reverse_backend_v2.sh
# or, if unavailable:
#   install_had_quadra_flat_accumulator_reverse_backend_v1.sh
#
# It intentionally keeps:
#   core/had/batch_directional_flat_accumulator.hpp
#   benchmarks/benchmark_batch_directional_flat_accumulator.cpp
#   benchmarks/benchmark_sparse_rw1_flat_accumulator_trace_check.cpp
#   run_batch_directional_flat_accumulator_benchmark.sh
#   run_sparse_rw1_flat_accumulator_trace_check.sh
#
# because those proved the flat accumulator architecture.

set -euo pipefail

target="core/had_quadra.hpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

backup=""
if ls .quadra_patch_backups/had_quadra.hpp.flat_accumulator_reverse_backend_v2.*.bak >/dev/null 2>&1; then
  backup="$(ls -t .quadra_patch_backups/had_quadra.hpp.flat_accumulator_reverse_backend_v2.*.bak | head -1)"
elif ls .quadra_patch_backups/had_quadra.hpp.flat_accumulator_reverse_backend.*.bak >/dev/null 2>&1; then
  backup="$(ls -t .quadra_patch_backups/had_quadra.hpp.flat_accumulator_reverse_backend.*.bak | head -1)"
else
  echo "ERROR: no flat-accumulator reverse backend backup found."
  echo "Available had_quadra backups:"
  ls -t .quadra_patch_backups/had_quadra.hpp.*.bak 2>/dev/null | head -20 || true
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.before_revert_flat_accumulator_reverse_backend.$(date +%Y%m%d_%H%M%S).bak"
cp "$backup" "$target"

echo
echo "Reverted core/had_quadra.hpp using:"
echo "  $backup"
echo
echo "Kept flat accumulator scaffold and trace-check benchmark files."
echo
echo "Run:"
echo "  ./run_had_quadra_nonzero_batch_directional_test.sh"
echo "  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10"
echo "  ./run_sparse_rw1_flat_accumulator_trace_check.sh 10"
echo
