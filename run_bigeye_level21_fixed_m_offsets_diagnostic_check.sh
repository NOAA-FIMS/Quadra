#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp"
OUT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs"
WORK="examples/NMFS/pifsc_bigeye_tuna/workflow"

mkdir -p "$WORK"

echo "== Confirm M offsets are fixed in objective =="
grep -n "log_m_young_offset = T(0.0)\|log_m_old_offset = T(0.0)\|fixed-M diagnostic" "$OBJ"

echo
echo "== Run Level 21 fixed-M-offset diagnostic =="
./run_bigeye_level21_age_based_m_check.sh

cp "$OUT/bigeye_level21_fit_summary.csv" "$WORK/bigeye_level21_fixed_m_offsets_fit_summary.csv"
cp "$OUT/bigeye_level21_gradient_by_parameter.csv" "$WORK/bigeye_level21_fixed_m_offsets_gradient_by_parameter.csv"
cp "$OUT/bigeye_level21_parameter_sanity_diagnostics.csv" "$WORK/bigeye_level21_fixed_m_offsets_parameter_sanity_diagnostics.csv"
cp "$OUT/bigeye_level21_age_comp_residual_diagnostics.txt" "$WORK/bigeye_level21_fixed_m_offsets_age_comp_residual_diagnostics.txt"

echo
echo "== Fixed-M fit summary =="
grep -E "objective|grad_norm|converged|log_m_young_offset|log_m_old_offset|init_number_multiplier_age_" "$OUT/bigeye_level21_fit_summary.csv"

echo
echo "== Fixed-M top gradients =="
grep -A25 -n "Top gradients" "$OUT/bigeye_level21_gradient_by_parameter.txt"

echo
echo "== Compare against current/free-M baseline if available =="
cat > "$WORK/bigeye_level21_fixed_m_offsets_interpretation.txt" <<'TXT'
Level 21 Fixed-M-Offsets Diagnostic
===================================

Purpose
-------
This diagnostic keeps the Level 21 parameter layout intact, but fixes
log_m_young_offset = 0 and log_m_old_offset = 0 inside the objective.

Interpretation rule
-------------------
If objective, residuals, and gradient norm improve relative to the free-M run,
then age-based M freedom was likely creating confounding or optimizer curvature
problems.

If objective worsens but gradient norm improves, M flexibility may be improving
fit but making optimization harder.

If both objective and gradient norm worsen, age-based M was probably helping,
and the remaining issue is more likely optimizer control, scaling, or another
model component.
TXT

echo "wrote: $WORK/bigeye_level21_fixed_m_offsets_interpretation.txt"
