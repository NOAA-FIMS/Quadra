#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"

ts="$(date +%Y%m%d_%H%M%S)"
cp "$LL" "$LL.before_bh_longline_weight_maturity_fix.$ts"
echo "backup: $LL.before_bh_longline_weight_maturity_fix.$ts"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/diagnostics/bigeye_longline_prediction_decomposition.hpp")
txt = p.read_text()

# The earlier patch inserted a phi0/spawning-biomass calculation before the
# local weight/maturity declarations existed. Add them immediately before the
# first spawning_biomass_from_numbers(n, weight, maturity) use if needed.
needle = "level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);"
idx = txt.find(needle)
if idx < 0:
    raise SystemExit("Could not find spawning_biomass_from_numbers(n, weight, maturity) in longline diagnostic")

prefix = txt[:idx]
block = ""
if "const auto weight = default_weight_at_age();" not in prefix:
    block += "    const auto weight = default_weight_at_age();\n"
if "const auto maturity = default_maturity_at_age();" not in prefix:
    block += "    const auto maturity = default_maturity_at_age();\n"

if block:
    # Insert at the start of the statement containing the needle.
    line_start = txt.rfind("\n", 0, idx) + 1
    txt = txt[:line_start] + block + txt[line_start:]
    print("inserted missing weight/maturity before first BH phi0 use")
else:
    print("weight/maturity already declared before first BH phi0 use")

p.write_text(txt)
PY

cat > run_bigeye_level23_longline_weight_maturity_fix_check.sh <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"

echo "== Longline diagnostic BH local context =="
grep -n "default_weight_at_age\|default_maturity_at_age\|spawning_biomass_from_numbers\|const double phi0\|expected_recruitment\|next\[0\]" "$LL"

echo
echo "== Build/run Level 23 BH diagnostic sync check =="
./run_bigeye_level23_bh_diag_compile_fix_check.sh
EOS
chmod +x run_bigeye_level23_longline_weight_maturity_fix_check.sh

echo
echo "Created: ./run_bigeye_level23_longline_weight_maturity_fix_check.sh"
echo "Run:"
echo "  ./run_bigeye_level23_longline_weight_maturity_fix_check.sh"
