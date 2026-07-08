#!/usr/bin/env bash
set -euo pipefail

REPORT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/reports/bigeye_fit_reports.hpp"

cp "$REPORT" "$REPORT.before_level21_fit_summary_loop_surgical_fix"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/reports/bigeye_fit_reports.hpp")
s = p.read_text()

lines = s.splitlines()
new_lines = []
for line in lines:
    if "init_log_number_dev_age_" in line or "init_number_multiplier_age_" in line:
        # Old/wrong: reporting selectivity offset as initial dev offset.
        line = re.sub(r'fit\.par\s*\[\s*kLonglineSelOffset\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', line)
        line = re.sub(r'par\s*\[\s*kLonglineSelOffset\s*\+\s*a\s*\]', 'par[3 + 2 + kAges + a]', line)
        line = re.sub(r'fit\.par\s*\[\s*level21_report_layout::kLonglineSelOffset\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', line)
        line = re.sub(r'par\s*\[\s*level21_report_layout::kLonglineSelOffset\s*\+\s*a\s*\]', 'par[3 + 2 + kAges + a]', line)
        line = re.sub(r'fit\.par\s*\[\s*kBaseFixed\s*\+\s*kAges\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', line)
        line = re.sub(r'par\s*\[\s*kBaseFixed\s*\+\s*kAges\s*\+\s*a\s*\]', 'par[3 + 2 + kAges + a]', line)
    new_lines.append(line)

s = "\n".join(new_lines) + ("\n" if s.endswith("\n") else "")

# Also handle the common pattern where a local init_dev variable is assigned once and then printed.
s = re.sub(
    r'(const\s+(?:double|auto)\s+init_dev\s*=\s*fit\.par\s*\[\s*)(?:kLonglineSelOffset|level21_report_layout::kLonglineSelOffset|kBaseFixed\s*\+\s*kAges)\s*(\+\s*a\s*\]\s*;)',
    r'\1(3 + 2 + kAges)\2',
    s
)
s = re.sub(
    r'(const\s+(?:double|auto)\s+init_dev\s*=\s*par\s*\[\s*)(?:kLonglineSelOffset|level21_report_layout::kLonglineSelOffset|kBaseFixed\s*\+\s*kAges)\s*(\+\s*a\s*\]\s*;)',
    r'\1(3 + 2 + kAges)\2',
    s
)

p.write_text(s)
print(f"patched {p}")
PY

echo "== Show init summary emit lines =="
grep -n "init_log_number_dev_age_\|init_number_multiplier_age_\|init_dev" "$REPORT" | head -100

echo
echo "Now run:"
echo "  ./run_bigeye_level21_fit_summary_layout_check.sh"
