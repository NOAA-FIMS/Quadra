#!/usr/bin/env bash
set -euo pipefail

LOG="${1:-level21_logdet_gradient_sign_probe.log}"

./run_bigeye_level21_age_based_m_check.sh > "$LOG" 2>&1

echo "wrote: $LOG"
echo
echo "== Profiled FD summary =="
grep -E "summary,(max_abs_analytic_minus_fd|max_abs_analytic_plus_fd|dot_analytic_fd|norm_analytic|norm_fd|cosine_analytic_fd|interpretation)" "$LOG" | tail -7 || true

echo
echo "== Gradient parts summary =="
grep -E "summary,(total_cosine|joint_cosine|logdet_cosine|total_norm_analytic|total_norm_fd|joint_norm_analytic|joint_norm_fd|logdet_norm_analytic|logdet_norm_fd)" "$LOG" | tail -9 || true

echo
echo "== Fit summary key rows =="
grep -E "^(objective|grad_norm|converged)," \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv 2>/dev/null || true
