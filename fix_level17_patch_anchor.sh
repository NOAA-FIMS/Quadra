#!/usr/bin/env bash
set -euo pipefail

PATCH="patch_bigeye_level17_juvenile_mortality_diagnostic.sh"

if [[ ! -f "$PATCH" ]]; then
  echo "ERROR: missing $PATCH in current directory."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$PATCH" "${PATCH}.before_anchor_fix.${STAMP}"

python3 - <<'PY'
from pathlib import Path

p = Path("patch_bigeye_level17_juvenile_mortality_diagnostic.sh")
s = p.read_text()

start = s.find("anchor = '    for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a)")
if start == -1:
    raise SystemExit("Could not find brittle anchor block start in patch script")

end_marker = 'driver.write_text(s)\n'
end = s.find(end_marker, start)
if end == -1:
    raise SystemExit("Could not find brittle anchor block end in patch script")
end += len(end_marker)

replacement = (
"insert = '    params.add({\"log_juvenile_m_multiplier\", 0.0,\\n"
"                quadra::ParameterTransform::Identity, false});\\n\\n'\n"
"if \"log_juvenile_m_multiplier\" not in s:\n"
"    candidates = [\n"
"        '    for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a) {\\\\n      params.add({\"init_log_number_dev_age_\"',\n"
"        '    for (int a = 0; a < kAges; ++a) {\\\\n      params.add({\"init_log_number_dev_age_\"',\n"
"        'init_log_number_dev_age_1',\n"
"    ]\n"
"    placed = False\n"
"    for anchor in candidates:\n"
"        idx = s.find(anchor)\n"
"        if idx != -1:\n"
"            line_start = s.rfind('\\n', 0, idx) + 1\n"
"            s = s[:line_start] + insert + s[line_start:]\n"
"            placed = True\n"
"            break\n"
"    if not placed:\n"
"        marker = 'logit_sel_purse_seine_age_'\n"
"        idx = s.find(marker)\n"
"        if idx != -1:\n"
"            block_start = s.rfind('    for (int a = 0;', 0, idx)\n"
"            if block_start != -1:\n"
"                s = s[:block_start] + insert + s[block_start:]\n"
"                placed = True\n"
"    if not placed:\n"
"        raise SystemExit('Could not find any parameter insertion point for log_juvenile_m_multiplier')\n"
"driver.write_text(s)\n"
)

s = s[:start] + replacement + s[end:]
p.write_text(s)
PY

echo "Fixed anchor logic in $PATCH"
echo
echo "Now rerun:"
echo "  ./patch_bigeye_level17_juvenile_mortality_diagnostic.sh"
echo "  ./run_bigeye_level17_juvenile_mortality_check.sh"
