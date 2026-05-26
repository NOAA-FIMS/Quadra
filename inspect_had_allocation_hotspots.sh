#!/usr/bin/env bash
set -euo pipefail

# inspect_had_allocation_hotspots.sh
#
# Locate exact had.h / Quadra allocation hot spots before patching internals.
# This script does not modify code.
#
# Run:
#   bash inspect_had_allocation_hotspots.sh
#
# Output:
#   quadra_had_allocation_hotspots.txt

OUT="${OUT:-quadra_had_allocation_hotspots.txt}"

{
echo "Quadra / had allocation hotspot inspection"
echo "Generated: $(date)"
echo

echo "== Compiler / platform =="
uname -a || true
echo

echo "== Candidate files containing had::NewAReal or NewAReal =="
grep -RIn "NewAReal" . \
  --exclude-dir=.git \
  --exclude-dir=.quadra_patch_backups \
  --exclude='*.trace*' \
  --exclude='*.xml' \
  2>/dev/null || true
echo

echo "== Candidate files containing ADGraph =="
grep -RIn "class ADGraph\|struct ADGraph\|ADGraph::\|g_ADGraph\|ADGraph" external core benchmarks \
  --exclude-dir=.git \
  --exclude-dir=.quadra_patch_backups \
  2>/dev/null || true
echo

echo "== Candidate files containing BTree::Insert or BTree =="
grep -RIn "BTree::Insert\|class BTree\|struct BTree\|BTree" external core benchmarks \
  --exclude-dir=.git \
  --exclude-dir=.quadra_patch_backups \
  2>/dev/null || true
echo

echo "== Direct heap allocation sites in external/had and core/autodiff/laplace =="
grep -RIn "new \|delete \|malloc\|free\|std::vector<\|push_back\|emplace_back" external/had core/autodiff core/laplace \
  --exclude-dir=.git \
  --exclude-dir=.quadra_patch_backups \
  2>/dev/null || true
echo

echo "== evaluate_random_effect_hessian definitions and call sites =="
grep -RIn "evaluate_random_effect_hessian" core benchmarks \
  --exclude-dir=.git \
  --exclude-dir=.quadra_patch_backups \
  2>/dev/null || true
echo

echo "== Laplace cached/evaluator files =="
find core benchmarks -type f \
  \( -iname '*laplace*' -o -iname '*hessian*' -o -iname '*autodiff*' \) \
  | sort
echo

echo "== Current workspace hooks =="
grep -RIn "HadQuadraWorkspace\|had_quadra_workspace" core benchmarks \
  --exclude-dir=.git \
  --exclude-dir=.quadra_patch_backups \
  2>/dev/null || true
echo

echo "== Make targets related to Laplace / workspace / had / arena =="
grep -n "benchmark-laplace\|workspace\|had\|arena" Makefile 2>/dev/null || true
echo

} > "${OUT}"

echo "Wrote ${OUT}"
echo
echo "Useful next command:"
echo "  sed -n '1,240p' ${OUT}"
echo
echo "Or upload/paste ${OUT} here."
