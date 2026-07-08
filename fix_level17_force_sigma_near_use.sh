#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic/objective/bigeye_quadra_objective.hpp"

if [[ ! -f "$OBJ" ]]; then
  echo "ERROR: missing $OBJ"
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$OBJ" "${OBJ}.before_force_sigma_near_use.${STAMP}"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic/objective/bigeye_quadra_objective.hpp")
s = p.read_text()

# Remove every existing declaration so there is exactly one.
s = re.sub(r"\n\s*const T sigma_log_juvenile_m_multiplier = T\(0\.50\);", "", s)

prior = "    nll = nll + T(0.5) * square_t(log_juvenile_m_multiplier / sigma_log_juvenile_m_multiplier);"
if prior not in s:
    raise SystemExit("Could not find juvenile M prior line")

# Put the declaration immediately before the prior. This is ugly but unambiguous.
s = s.replace(
    prior,
    "    const T sigma_log_juvenile_m_multiplier = T(0.50);\n" + prior,
    1,
)

p.write_text(s)
PY

echo "Forced sigma_log_juvenile_m_multiplier declaration immediately before use."
echo
echo "Run:"
echo "  grep -n \"sigma_log_juvenile\\|log_juvenile_m_multiplier\" \"$OBJ\""
echo "  ./run_bigeye_level17_juvenile_mortality_check.sh"
