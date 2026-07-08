#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/objective/bigeye_quadra_objective.hpp"

echo "== Confirm Level 23 BH uses computed unfished phi0 =="
grep -n "unfished_spawning_biomass\|const T phi0\|expected_recruitment\|next\[0\].*expected_recruitment" "${OBJ}" || true

echo
echo "== Confirm placeholder phi0 is gone =="
if grep -n "phi0 = T(1.0)" "${OBJ}"; then
  echo "ERROR: placeholder phi0 remains" >&2
  exit 1
fi

echo
echo "== Run Level 23 computed-phi0 BH check =="
if [[ -x ./run_bigeye_level23_longline_selectivity_smoothness_check.sh ]]; then
  ./run_bigeye_level23_longline_selectivity_smoothness_check.sh
elif [[ -x ./run_bigeye_level23_longline_selectivity_smoothness_scan.sh ]]; then
  ./run_bigeye_level23_longline_selectivity_smoothness_scan.sh
else
  echo "No known Level 23 runner found."
  find examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan -maxdepth 3 -type f -name "run*.sh" -print | sort
  exit 1
fi

echo
echo "== Level 21 vs Level 23 computed-phi0 BH comparison =="
grep -E "objective|grad_norm|converged" \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv \
  examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs/bigeye_level23_fit_summary.csv || true

echo
echo "== Best Level 23 scan rows if summary exists =="
summary="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs/bigeye_level23_smoothness_scan_summary.csv"
if [[ -f "${summary}" ]]; then
  python3 - <<'PY'
from pathlib import Path
import csv
p = Path("examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs/bigeye_level23_smoothness_scan_summary.csv")
rows = list(csv.DictReader(p.open()))
rows = sorted(rows, key=lambda r: float(r["objective"]))
for r in rows[:7]:
    print(",".join([r.get(k,"") for k in ["ll_sigma","smooth_lambda","objective","grad_norm","converged","longline_mean_abs","purse_seine_mean_abs","longline_selectivity_prior_nll","initial_numbers_prior_nll","recruitment_prior_nll"]]))
PY
else
  echo "No summary CSV found at ${summary}"
fi
