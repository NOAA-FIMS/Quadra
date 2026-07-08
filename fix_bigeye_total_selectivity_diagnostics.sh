#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna"

python3 - <<'PY'
from pathlib import Path
import re

ROOT = Path("examples/NMFS/pifsc_bigeye_tuna")
levels = [
    "level20_longline_selectivity_regularization_scan",
    "level21_age_based_natural_mortality_diagnostic",
    "level23_longline_selectivity_smoothness_scan",
]

changed = []

def backup(p: Path):
    b = p.with_name(p.name + ".before_total_selectivity_diag_fix")
    if not b.exists():
        b.write_text(p.read_text())

for level in levels:
    diag = ROOT / level / "diagnostics"
    if not diag.exists():
        continue
    for p in diag.glob("*.hpp"):
        s0 = p.read_text()
        s = s0

        # Fix calculation blocks that still use avg selectivity for total F/Z.
        s = s.replace(
            "avg_sel[i] = 0.5 * (sel_longline[i] + sel_purse_seine[i]);",
            "avg_sel[i] = sel_longline[i] + sel_purse_seine[i];"
        )
        s = s.replace(
            "avg_sel[i] = T(0.5) * (sel_longline[i] + sel_purse_seine[i]);",
            "avg_sel[i] = sel_longline[i] + sel_purse_seine[i];"
        )

        # In plus-group audit, previous patch renamed the local variable but not print statements.
        s = s.replace(" << avg_sel << ", " << total_sel << ")
        s = s.replace("<< avg_sel <<", "<< total_sel <<")
        s = s.replace("avg_sel,z,annual_survival", "total_sel,z,annual_survival")
        s = s.replace("avg_sel,z,annual_survival", "total_sel,z,annual_survival")

        # If comments/text label avg_sel, make it clear this is total fleet selectivity.
        s = s.replace("avg_sel", "total_sel")

        if s != s0:
            backup(p)
            p.write_text(s)
            changed.append(str(p))

print("Patched diagnostic/report total-selectivity naming and calculations:")
for p in changed:
    print("  " + p)
if not changed:
    print("  none")
PY

cat > run_bigeye_total_selectivity_f_fix_check_v2.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna"

echo "== Check for remaining fleet-average F patterns =="
grep -R "0\.5.*sel_longline.*sel_purse_seine\|T(0\.5).*sel_longline.*sel_purse_seine" -n \
  "$ROOT"/level20_longline_selectivity_regularization_scan \
  "$ROOT"/level21_age_based_natural_mortality_diagnostic \
  "$ROOT"/level23_longline_selectivity_smoothness_scan \
  --include='*.hpp' --include='*.cpp' || true

echo
echo "== Check total-selectivity F patterns =="
grep -R "total_sel.*sel_longline.*sel_purse_seine\|total_sel_prev.*sel_longline.*sel_purse_seine\|total_sel_last.*sel_longline.*sel_purse_seine" -n \
  "$ROOT"/level20_longline_selectivity_regularization_scan \
  "$ROOT"/level21_age_based_natural_mortality_diagnostic \
  "$ROOT"/level23_longline_selectivity_smoothness_scan \
  --include='*.hpp' --include='*.cpp' | head -160

echo
echo "== Rerun Level 21 age-M check =="
./run_bigeye_level21_age_based_m_check.sh
SH

chmod +x run_bigeye_total_selectivity_f_fix_check_v2.sh

echo
echo "Run:"
echo "  ./run_bigeye_total_selectivity_f_fix_check_v2.sh"
