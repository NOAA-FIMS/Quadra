#!/usr/bin/env bash
set -euo pipefail

L17="examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic"

if [[ ! -d "$L17" ]]; then
  echo "ERROR: missing $L17."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
for f in \
  "$L17/objective/bigeye_quadra_objective.hpp" \
  "$L17/diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
do
  [[ -f "$f" ]] && cp "$f" "${f}.before_scope_fix.${STAMP}"
done

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic")

# Objective: guarantee sigma is declared before juvenile prior.
p = root / "objective/bigeye_quadra_objective.hpp"
s = p.read_text()

# Remove all existing juvenile sigma declarations so we control placement.
s = re.sub(r"\n\s*const T sigma_log_juvenile_m_multiplier = T\(0\.50\);", "", s)

# Insert immediately after sigma_ps_sel_dev declaration.
needle = "const T sigma_ps_sel_dev = T(1.0);"
if needle not in s:
    raise SystemExit("Could not find sigma_ps_sel_dev in objective")
s = s.replace(
    needle,
    needle + "\n    const T sigma_log_juvenile_m_multiplier = T(0.50);",
    1,
)

p.write_text(s)

# Age-comp residual diagnostics: guarantee kJuvenileMIndex is declared before first use.
p = root / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
s = p.read_text()

# Normalize layout constants.
s = s.replace("constexpr int kBaseFixed = 6;", "constexpr int kBaseFixed = 5;")

# Remove any duplicated kJuvenileMIndex lines.
s = re.sub(r"\n\s*constexpr int kJuvenileMIndex = kBaseFixed \+ kInitialDevs;", "", s)

# Insert directly after kInitialDevs.
needle = "constexpr int kInitialDevs = kAges;"
if needle not in s:
    raise SystemExit("Could not find kInitialDevs in age-comp diagnostics")
s = s.replace(
    needle,
    needle + "\n  constexpr int kJuvenileMIndex = kBaseFixed + kInitialDevs;",
    1,
)

# Ensure PS offset uses juvenile index.
s = re.sub(
    r"constexpr int kPurseSeineSelOffset = [^;]+;",
    "constexpr int kPurseSeineSelOffset = kJuvenileMIndex + 1;",
    s,
    count=1,
)

p.write_text(s)
PY

echo "Fixed Level 17 scope/declaration order."
echo
echo "Run:"
echo "  ./run_bigeye_level17_juvenile_mortality_check.sh"
