#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Find all writers of bigeye_level21_fit_summary.csv and init rows =="
grep -R "bigeye_level21_fit_summary\|init_number_multiplier_age_" -n "$L21" \
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
    if "init_number_multiplier_age_" not in s:
        continue

    orig = s
    backup = p.with_name(p.name + ".before_level21_init_multiplier_emit_override")
    if not backup.exists():
        backup.write_text(s)

    # Replace stream lines that emit fit-summary init multipliers.
    # Level 21 true init dev offset = 3 base + 2 M params + kAges LL selectivity params.
    # Do not touch the parameter sanity table because it is already correct.
    lines = s.splitlines()
    out = []
    for line in lines:
        if "init_number_multiplier_age_" in line and "<<" in line:
            indent = line[:len(line) - len(line.lstrip())]
            stream = "csv"
            m = re.match(r'\s*([A-Za-z_][A-Za-z0-9_]*)\s*<<', line)
            if m:
                stream = m.group(1)
            vector_name = "fit.par" if "fit.par" in s else "par"
            line = (
                indent
                + f'{stream} << "init_number_multiplier_age_" << (a + 1) << "," '
                + f'<< std::exp({vector_name}[3 + 2 + kAges + a]) << "\\n";'
            )
        out.append(line)
    s = "\n".join(out) + ("\n" if orig.endswith("\n") else "")

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
echo "== Emit sites after override =="
grep -R "init_number_multiplier_age_" -n "$L21" \
  --include='*.hpp' --include='*.cpp' \
  | grep -v '\.before_' || true

cat > run_bigeye_level21_init_multiplier_emit_override_check.sh <<'SH'
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
echo "Expected corrected fit summary:"
echo "  age1 ~1.029, age2 ~0.999, age8 ~0.724, age9 ~0.650, age10 ~0.671"
SH

chmod +x run_bigeye_level21_init_multiplier_emit_override_check.sh

echo
echo "Run:"
echo "  ./run_bigeye_level21_init_multiplier_emit_override_check.sh"
