#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
DIAG_DIR="$LVL/diagnostics"
ts="$(date +%Y%m%d_%H%M%S)"

files=(
  "$DIAG_DIR/bigeye_age_comp_residual_diagnostics.hpp"
  "$DIAG_DIR/bigeye_purse_seine_prediction_decomposition.hpp"
  "$DIAG_DIR/bigeye_longline_prediction_decomposition.hpp"
  "$DIAG_DIR/bigeye_level21_plus_group_audit.hpp"
)

for f in "${files[@]}"; do
  [[ -f "$f" ]] || continue
  cp "$f" "$f.before_level21_phi0_diag_fix.$ts"
  echo "backup: $f.before_level21_phi0_diag_fix.$ts"
done

python3 - <<'PY'
from pathlib import Path
import re

lvl = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")
diag_dir = lvl / "diagnostics"

files = [
    diag_dir / "bigeye_age_comp_residual_diagnostics.hpp",
    diag_dir / "bigeye_purse_seine_prediction_decomposition.hpp",
    diag_dir / "bigeye_longline_prediction_decomposition.hpp",
    diag_dir / "bigeye_level21_plus_group_audit.hpp",
]

phi0_block = """    const double phi0 =
        level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity) /
        std::max(r0, level21_bh_sync::bh_min_positive());
"""

for p in files:
    if not p.exists():
        continue

    s = p.read_text()
    if "level21_bh_sync::beverton_holt_recruitment" not in s:
        continue
    if "const double phi0 =" in s:
        print(f"phi0 already present: {p}")
        continue

    target = re.search(
        r'(\s*const double spawning_biomass\s*=\s*\n\s*level21_bh_sync::spawning_biomass_from_numbers\(n,\s*weight,\s*maturity\);)',
        s
    )
    if target:
        s = s[:target.start()] + "\n" + phi0_block + s[target.start():]
        p.write_text(s)
        print(f"inserted phi0 before BH recruitment block: {p}")
    else:
        raise SystemExit(f"ERROR: could not find BH spawning_biomass block in {p}")
PY

cat > run_bigeye_level21_phi0_diag_fix_check.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Level 21 diagnostic phi0 context =="
grep -R "const double phi0\|beverton_holt_recruitment\|next\[0\]" -n \
  "$LVL/diagnostics" --include='*.hpp' | grep -v '\.before_' || true

echo
echo "== Build/run Level 21 BH check =="
./run_bigeye_level21_age_based_m_check.sh
EOF

chmod +x run_bigeye_level21_phi0_diag_fix_check.sh

echo
echo "Created:"
echo "  ./run_bigeye_level21_phi0_diag_fix_check.sh"
echo
echo "Run:"
echo "  ./run_bigeye_level21_phi0_diag_fix_check.sh"
