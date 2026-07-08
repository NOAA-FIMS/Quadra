#!/usr/bin/env bash
set -euo pipefail

REPORT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/reports/bigeye_fit_reports.hpp"

if [[ ! -f "$REPORT" ]]; then
  echo "ERROR: report file not found: $REPORT"
  exit 1
fi

cp "$REPORT" "$REPORT.before_level21_report_constexpr_alias_cleanup"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/reports/bigeye_fit_reports.hpp")
s = p.read_text()

# Replace qualified constexpr aliases with literal expressions so they are true constant expressions.
repls = {
    r"constexpr int kLonglineSelOffset\s*=\s*level21_report_layout::kLonglineSelOffset\s*;":
        "constexpr int kLonglineSelOffset = kBaseFixed + 2;",
    r"constexpr int kInitialDevOffset\s*=\s*level21_report_layout::kInitialDevOffset\s*;":
        "constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;",
    r"constexpr int kPurseSeineSelOffset\s*=\s*level21_report_layout::kPurseSeineSelOffset\s*;":
        "constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;",
    r"constexpr int kRecruitmentOffset\s*=\s*level21_report_layout::kRecruitmentOffset\s*;":
        "constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;",
}
for pat, repl in repls.items():
    s = re.sub(pat, repl, s)

# If an inserted namespace-level layout block exists and is now redundant, leave it alone;
# it is harmless. The key is removing the qualified constexpr aliases causing compile failure.

# Prevent accidental old Level 20 offset if present.
s = re.sub(r"constexpr int kLonglineSelOffset\s*=\s*kBaseFixed\s*;",
           "constexpr int kLonglineSelOffset = kBaseFixed + 2;", s)

p.write_text(s)
print(f"cleaned constexpr aliases in {p}")
PY

echo
echo "== Check bad qualified constexpr aliases are gone =="
grep -n "constexpr int k.*= level21_report_layout::" "$REPORT" || true

echo
echo "== Relevant layout constants =="
grep -n "kBaseFixed\|kLonglineSelOffset\|kInitialDevOffset\|kPurseSeineSelOffset\|kRecruitmentOffset" "$REPORT" | head -40

echo
echo "Now run:"
echo "  ./run_bigeye_level21_fit_summary_layout_check.sh"
