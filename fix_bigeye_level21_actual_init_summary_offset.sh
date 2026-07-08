#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Locate actual fit-summary init multiplier emit site =="
grep -R "init_number_multiplier_age_" -n "$L21" \
  --include='*.hpp' --include='*.cpp' \
  | grep -v '\.before_' || true

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")
files = list(root.rglob("*.hpp")) + list(root.rglob("*.cpp"))

changed = []
for p in files:
    if ".before_" in p.name:
        continue
    s = p.read_text(errors="replace")
    if "init_number_multiplier_age_" not in s and "init_log_number_dev_age_" not in s:
        continue

    orig = s
    backup = p.with_name(p.name + ".before_actual_init_summary_offset_fix")
    if not backup.exists():
        backup.write_text(s)

    # In any file that emits init-number rows, force Level 21 init offset:
    # 3 base fixed + 2 M params + kAges longline selectivity devs.
    s = re.sub(
        r'(const\s+(?:double|auto)\s+init_dev\s*=\s*)(?:fit\.)?par\s*\[\s*(?:kLonglineSelOffset|level21_report_layout::kLonglineSelOffset|kBaseFixed\s*\+\s*kAges)\s*\+\s*a\s*\]\s*;',
        r'\1fit.par[3 + 2 + kAges + a];',
        s
    )

    s = re.sub(r'fit\.par\s*\[\s*kLonglineSelOffset\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', s)
    s = re.sub(r'par\s*\[\s*kLonglineSelOffset\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', s)
    s = re.sub(r'fit\.par\s*\[\s*level21_report_layout::kLonglineSelOffset\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', s)
    s = re.sub(r'par\s*\[\s*level21_report_layout::kLonglineSelOffset\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', s)
    s = re.sub(r'fit\.par\s*\[\s*kBaseFixed\s*\+\s*kAges\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', s)
    s = re.sub(r'par\s*\[\s*kBaseFixed\s*\+\s*kAges\s*\+\s*a\s*\]', 'fit.par[3 + 2 + kAges + a]', s)

    if s != orig:
        p.write_text(s)
        changed.append(str(p))

print("changed files:")
for c in changed:
    print("  " + c)
if not changed:
    print("  none")
PY

echo
echo "== Re-show init emit sites after patch =="
grep -R "init_number_multiplier_age_\|init_log_number_dev_age_\|init_dev" -n "$L21" \
  --include='*.hpp' --include='*.cpp' \
  | grep -v '\.before_' | head -120 || true

cat > run_bigeye_level21_actual_init_summary_offset_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

./run_bigeye_level21_age_based_m_check.sh

FIT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv"
SAN="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_parameter_sanity_diagnostics.txt"

echo
echo "== Fit summary init multipliers =="
grep "init_number_multiplier_age_" "$FIT"

echo
echo "== Sanity table init multipliers =="
grep -A12 -n "Age-specific parameters" "$SAN"

echo
echo "Expected age 8/9/10:"
echo "  age8 ~0.724, age9 ~0.650, age10 ~0.671"
SH

chmod +x run_bigeye_level21_actual_init_summary_offset_check.sh

echo
echo "Run:"
echo "  ./run_bigeye_level21_actual_init_summary_offset_check.sh"
