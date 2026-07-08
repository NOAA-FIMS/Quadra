#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Objective layout constants =="
grep -n "kBaseFixed\|kMParamOffset\|kMParams\|kLonglineSelOffset\|kInitialDevOffset\|kPurseSeineSelOffset\|kRecruitmentOffset" \
  "$L21/objective/bigeye_quadra_objective.hpp" | head -80

echo
echo "== Report/diagnostic stale hard-coded offset scan =="
grep -R "fit\.par\[[0-9][0-9]* *+ *a\]\|fit\.par\[[0-9][0-9]* *+ *static_cast" -n \
  "$L21/reports" "$L21/diagnostics" || true

echo
echo "== M-at-age usage scan =="
grep -R "m_at_age\|log_m_young_offset\|log_m_old_offset\|m_young\|m_old" -n \
  "$L21/objective" "$L21/reports" "$L21/diagnostics" | head -160 || true
