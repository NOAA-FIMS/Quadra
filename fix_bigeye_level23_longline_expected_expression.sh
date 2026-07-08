#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"

ts="$(date +%Y%m%d_%H%M%S)"
cp "$LL" "$LL.before_bh_longline_expected_expression_fix.$ts"
echo "backup: $LL.before_bh_longline_expected_expression_fix.$ts"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/diagnostics/bigeye_longline_prediction_decomposition.hpp")
txt = p.read_text()

# Remove the bad declarations that were inserted in the middle of an expression.
txt = txt.replace(
    "    const auto weight = default_weight_at_age();\n"
    "    const auto maturity = default_maturity_at_age();\n"
    "      level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);",
    "      level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);"
)

# Insert declarations immediately after r0 is declared, where they are valid statements
# and in scope for phi0 and the yearly loop.
m = re.search(r'^(\s*const double r0\s*=.*;\n)', txt, flags=re.M)
if not m:
    raise SystemExit("Could not find const double r0 declaration in longline diagnostic")

insert_at = m.end()
prefix = txt[:insert_at]
suffix = txt[insert_at:]

decl = ""
if "const auto weight = default_weight_at_age();" not in suffix[:500]:
    decl += "  const auto weight = default_weight_at_age();\n"
if "const auto maturity = default_maturity_at_age();" not in suffix[:500]:
    decl += "  const auto maturity = default_maturity_at_age();\n"

if decl:
    txt = prefix + decl + suffix
    print("inserted weight/maturity after r0 declaration")
else:
    print("weight/maturity already present after r0 declaration")

p.write_text(txt)
PY

cat > run_bigeye_level23_longline_expected_expression_fix_check.sh <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"

echo "== Longline BH context =="
grep -n "const double r0\|default_weight_at_age\|default_maturity_at_age\|const double phi0\|spawning_biomass_from_numbers\|expected_recruitment\|next\[0\]" "$LL"

echo
echo "== Build/run Level 23 BH check =="
./run_bigeye_level23_bh_diag_compile_fix_check.sh
EOS
chmod +x run_bigeye_level23_longline_expected_expression_fix_check.sh

echo
echo "Created: ./run_bigeye_level23_longline_expected_expression_fix_check.sh"
echo "Run:"
echo "  ./run_bigeye_level23_longline_expected_expression_fix_check.sh"
