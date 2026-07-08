#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
DIAG="$L21/diagnostics"

if [[ ! -d "$DIAG" ]]; then
  echo "ERROR: missing $DIAG"
  exit 1
fi

python3 - <<'PY'
from pathlib import Path
import re

diag = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/diagnostics")

targets = [
    diag / "bigeye_age_comp_residual_diagnostics.hpp",
    diag / "bigeye_purse_seine_prediction_decomposition.hpp",
    diag / "bigeye_longline_prediction_decomposition.hpp",
    diag / "bigeye_initial_numbers_diagnostics.hpp",
    diag / "bigeye_level20_parameter_sanity_diagnostics.hpp",
]

def backup(p: Path):
    if p.exists():
        b = p.with_name(p.name + ".before_level21_layout_qualification_cleanup")
        if not b.exists():
            b.write_text(p.read_text())

# Lines like this are illegal inside functions:
# constexpr int level21_m_helpers::kLonglineSelOffset = ...
bad_local_def = re.compile(
    r'^\s*constexpr\s+int\s+level21_m_helpers::k[A-Za-z0-9_]+\s*=.*?;\s*\n',
    re.MULTILINE,
)

# Also remove any remaining unqualified local layout constexpr definitions from previous patch attempts.
bad_unqualified_def = re.compile(
    r'^\s*constexpr\s+int\s+k(?:BaseFixed|MParamOffset|MParams|LonglineSelOffset|LonglineSelDevs|InitialDevOffset|InitialDevs|PurseSeineSelOffset|PurseSeineSelDevs|RecruitmentOffset)\s*=.*?;\s*\n',
    re.MULTILINE,
)

for p in targets:
    if not p.exists():
        continue
    backup(p)
    s = p.read_text()
    s = bad_local_def.sub('', s)
    s = bad_unqualified_def.sub('', s)

    # Repair accidental double qualification.
    s = s.replace('level21_m_helpers::level21_m_helpers::', 'level21_m_helpers::')

    # Make sure helper include is present exactly once.
    inc = '#include "bigeye_level21_m_at_age_helpers.hpp"\n'
    s = s.replace(inc, '')
    if '#pragma once\n' in s:
        s = s.replace('#pragma once\n', '#pragma once\n' + inc, 1)
    else:
        s = inc + s

    p.write_text(s)

print("Cleaned illegal qualified local constexpr layout definitions.")
print("Now rerun:")
print("  ./run_bigeye_level21_age_based_m_check.sh")
PY
