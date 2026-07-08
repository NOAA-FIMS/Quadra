#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
AGE="$ROOT/quadra/bigeye_age_structured.hpp"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"
PS="$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp"
AC="$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp"

ts="$(date +%Y%m%d_%H%M%S)"
for f in "$LL" "$PS" "$AC"; do
  cp "$f" "$f.before_bh_diag_missing_locals_fix.$ts"
  echo "backup: $f.before_bh_diag_missing_locals_fix.$ts"
done

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan")
diags = [
    root / "diagnostics/bigeye_longline_prediction_decomposition.hpp",
    root / "diagnostics/bigeye_purse_seine_prediction_decomposition.hpp",
    root / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp",
]

def ensure_before_first(txt, marker, block):
    if block.strip() in txt:
        return txt, 0
    idx = txt.find(marker)
    if idx < 0:
        return txt, -1
    return txt[:idx] + block + txt[idx:], 1

for p in diags:
    txt = p.read_text()

    # Do not duplicate central helper. Diagnostics should use the helper from
    # bigeye_age_structured.hpp, which is already included through the objective path.
    txt = re.sub(
        r'\n#ifndef LEVEL23_BH_SYNC_HELPERS_DEFINED\n#define LEVEL23_BH_SYNC_HELPERS_DEFINED\nnamespace level23_bh_sync \{.*?\n#endif // LEVEL23_BH_SYNC_HELPERS_DEFINED\n',
        "\n", txt, flags=re.S
    )
    txt = re.sub(
        r'\nnamespace level23_bh_sync \{.*?\n\} // namespace level23_bh_sync\n',
        "\n", txt, flags=re.S
    )

    # If the BH diagnostic block references weight/maturity but this diagnostic
    # never declared them, declare them immediately before the first BH use.
    if "spawning_biomass_from_numbers(n, weight, maturity)" in txt:
        insert_marker = "const double spawning_biomass ="
        if "const auto weight = default_weight_at_age();" not in txt:
            txt, n = ensure_before_first(
                txt,
                insert_marker,
                "    const auto weight = default_weight_at_age();\n"
            )
            print(("inserted" if n == 1 else "WARNING: could not insert") + f" weight in {p}")
        if "const auto maturity = default_maturity_at_age();" not in txt:
            txt, n = ensure_before_first(
                txt,
                insert_marker,
                "    const auto maturity = default_maturity_at_age();\n"
            )
            print(("inserted" if n == 1 else "WARNING: could not insert") + f" maturity in {p}")

    # If phi0 is referenced but missing, compute a local per-current-state proxy
    # directly before expected_recruitment. This is a compile/diagnostic sync fix;
    # the objective remains the source of truth for the actual likelihood.
    if "beverton_holt_recruitment(spawning_biomass, r0, steepness, phi0)" in txt and "const double phi0 =" not in txt:
        marker = "const double expected_recruitment ="
        block = (
            "    const double phi0 =\n"
            "        level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity) /\n"
            "        std::max(r0, level23_bh_sync::bh_min_positive());\n"
        )
        txt, n = ensure_before_first(txt, marker, block)
        print(("inserted" if n == 1 else "WARNING: could not insert") + f" phi0 in {p}")

    p.write_text(txt)
PY

# Patch the checker: count namespace openings only, not the closing comment.
python3 - <<'PY'
from pathlib import Path
p = Path("run_bigeye_level23_bh_diag_compile_fix_check.sh")
if p.exists():
    txt = p.read_text()
    txt = txt.replace('grep -R "namespace level23_bh_sync" -n \\', 'grep -R "namespace level23_bh_sync {" -n \\')
    p.write_text(txt)
    print("patched checker namespace count")
PY

cat > run_bigeye_level23_bh_diag_missing_locals_fix_check.sh <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"

echo "== BH helper namespace openings =="
grep -R "namespace level23_bh_sync {" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

n=$(grep -R "namespace level23_bh_sync {" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' | wc -l | tr -d ' ')

if [[ "$n" != "1" ]]; then
  echo "ERROR: expected exactly one BH helper namespace opening, found $n" >&2
  exit 1
fi

echo
echo "== Remaining BH local-symbol references =="
grep -R "const auto weight = default_weight_at_age\|const auto maturity = default_maturity_at_age\|const double phi0\|expected_recruitment\|next\[0\]" -n \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

echo
echo "== Run previous full checker =="
./run_bigeye_level23_bh_diag_compile_fix_check.sh
EOS
chmod +x run_bigeye_level23_bh_diag_missing_locals_fix_check.sh

echo
echo "Created: ./run_bigeye_level23_bh_diag_missing_locals_fix_check.sh"
echo "Run:"
echo "  ./run_bigeye_level23_bh_diag_missing_locals_fix_check.sh"
