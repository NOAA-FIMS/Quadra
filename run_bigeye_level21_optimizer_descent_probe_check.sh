#!/usr/bin/env bash
set -euo pipefail

echo "== Confirm optimizer descent probe was inserted =="
grep -n "Optimizer Descent Probe\|fx_minus_alpha_g\|fx_plus_alpha_g\|LBFGSObjective" \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp || true

echo
echo "== Build/run Level 21 with optimizer descent probe =="
./run_bigeye_level21_age_based_m_check.sh | tee examples/NMFS/pifsc_bigeye_tuna/workflow/bigeye_level21_optimizer_descent_probe_run.log

echo
echo "== Extract descent probe output =="
awk '
  /Bigeye Level 21 Optimizer Descent Probe/ {show=1}
  show {print}
  show && /Interpretation:/ {seen_interp=1}
  seen_interp && /^$/ {exit}
' examples/NMFS/pifsc_bigeye_tuna/workflow/bigeye_level21_optimizer_descent_probe_run.log
