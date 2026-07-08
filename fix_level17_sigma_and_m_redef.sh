#!/usr/bin/env bash
set -euo pipefail

L17="examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic"

if [[ ! -d "$L17" ]]; then
  echo "ERROR: missing $L17."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$L17/objective/bigeye_quadra_objective.hpp" \
   "$L17/objective/bigeye_quadra_objective.hpp.before_sigma_fix.${STAMP}"
cp "$L17/diagnostics/bigeye_initial_numbers_diagnostics.hpp" \
   "$L17/diagnostics/bigeye_initial_numbers_diagnostics.hpp.before_m_redef_fix.${STAMP}"

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic")

# Fix missing sigma_log_juvenile_m_multiplier declaration.
p = root / "objective/bigeye_quadra_objective.hpp"
s = p.read_text()

if "const T sigma_log_juvenile_m_multiplier" not in s:
    # Place it immediately after sigma_ps_sel_dev if available.
    s = s.replace(
        "const T sigma_ps_sel_dev = T(1.0);",
        "const T sigma_ps_sel_dev = T(1.0);\n    const T sigma_log_juvenile_m_multiplier = T(0.50);",
        1,
    )

# If prior line appears before sigma declarations because of earlier patching,
# move it to after the base priors and after sigma exists.
prior = "    nll = nll + T(0.5) * square_t(log_juvenile_m_multiplier / sigma_log_juvenile_m_multiplier);\n"
if prior in s:
    s = s.replace(prior, "")
    anchor = "    nll = nll + normal_prior(log_sel_slope_longline, std::log(1.2), 0.35);\n"
    if anchor in s:
        s = s.replace(anchor, anchor + prior, 1)
    else:
        raise SystemExit("Could not find slope prior anchor for juvenile M prior")
else:
    anchor = "    nll = nll + normal_prior(log_sel_slope_longline, std::log(1.2), 0.35);\n"
    s = s.replace(anchor, anchor + prior, 1)

p.write_text(s)

# Remove stale duplicate m definition in initial numbers diagnostics.
p = root / "diagnostics/bigeye_initial_numbers_diagnostics.hpp"
s = p.read_text()

# Remove log_m/m lines if adult_m/juvenile_m/m already exist.
if "const double adult_m = 0.45;" in s:
    s = re.sub(r"\n\s*const double log_m = std::log\(0\.45\);", "", s)
    s = re.sub(r"\n\s*const double m = std::exp\(log_m\);", "", s)

p.write_text(s)
PY

echo "Fixed Level 17 sigma declaration and duplicate m definition."
echo
echo "Run:"
echo "  ./run_bigeye_level17_juvenile_mortality_check.sh"
