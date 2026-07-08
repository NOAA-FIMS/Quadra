#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
AGE="$ROOT/quadra/bigeye_age_structured.hpp"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"
PS="$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp"
AC="$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp"

ts="$(date +%Y%m%d_%H%M%S)"

for f in "$AGE" "$LL" "$PS" "$AC"; do
  [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
  cp "$f" "$f.before_bh_diag_compile_fix.$ts"
  echo "backup: $f.before_bh_diag_compile_fix.$ts"
done

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan")
age = root / "quadra/bigeye_age_structured.hpp"
diags = [
    root / "diagnostics/bigeye_longline_prediction_decomposition.hpp",
    root / "diagnostics/bigeye_purse_seine_prediction_decomposition.hpp",
    root / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp",
]

# The helper belongs in bigeye_age_structured.hpp only. Diagnostics include that
# indirectly, so duplicate helper namespaces in diagnostics cause same-TU redefinition.
helper_pat = re.compile(
    r'\nnamespace level23_bh_sync \{\n.*?\n\} // namespace level23_bh_sync\n',
    re.S
)

for p in diags:
    txt = p.read_text()
    txt2, n = helper_pat.subn("\n", txt)
    if n:
        print(f"removed duplicate BH helper from {p}: {n}")
    txt = txt2

    # Ensure local static inputs exist in the diagnostic/report function.
    # Insert after weight if present.
    if "const auto maturity = default_maturity_at_age();" not in txt:
        txt, n = re.subn(
            r'(const auto weight = default_weight_at_age\(\);\n)',
            r'\1  const auto maturity = default_maturity_at_age();\n',
            txt,
            count=1,
        )
        if n:
            print(f"inserted maturity in {p}")
        else:
            print(f"WARNING: no weight declaration found in {p}")

    # Ensure steepness exists after r0.
    if "const double steepness = 0.75;" not in txt:
        txt, n = re.subn(
            r'(const double r0 = std::exp\([^;]+;\n)',
            r'\1  const double steepness = 0.75;\n',
            txt,
            count=1,
        )
        if n:
            print(f"inserted steepness in {p}")
        else:
            print(f"WARNING: no r0 declaration found in {p}")

    # Ensure phi0 exists after n initialization. Some diagnostics use explicit loops
    # rather than unfished_equilibrium_numbers; try both placements.
    if "const double phi0 =" not in txt:
        txt, n = re.subn(
            r'(std::array<double, kAges> n = unfished_equilibrium_numbers\([^;]+;\n)',
            r'\1'
            r'  const double unfished_spawning_biomass =\n'
            r'      level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);\n'
            r'  const double phi0 =\n'
            r'      unfished_spawning_biomass / std::max(r0, level23_bh_sync::bh_min_positive());\n',
            txt,
            count=1,
        )
        if not n:
            # fallback: after initial n array declaration/initialization block before the year loop
            txt, n = re.subn(
                r'(\n\s*for \(std::size_t t = 0; t < years\.size\(\); \+\+t\) \{\n)',
                r'\n  const double unfished_spawning_biomass =\n'
                r'      level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);\n'
                r'  const double phi0 =\n'
                r'      unfished_spawning_biomass / std::max(r0, level23_bh_sync::bh_min_positive());\n'
                r'\1',
                txt,
                count=1,
            )
        if n:
            print(f"inserted phi0 in {p}")
        else:
            print(f"WARNING: no phi0 insertion point found in {p}")

    p.write_text(txt)

# Make the central helper include-guarded inside the namespace too, so even if a
# diagnostic later reintroduces it, this symbol block is harmlessly skipped.
txt = age.read_text()
if "LEVEL23_BH_SYNC_HELPERS_DEFINED" not in txt and "namespace level23_bh_sync" in txt:
    txt = txt.replace(
        "namespace level23_bh_sync {\n",
        "#ifndef LEVEL23_BH_SYNC_HELPERS_DEFINED\n#define LEVEL23_BH_SYNC_HELPERS_DEFINED\nnamespace level23_bh_sync {\n",
        1
    )
    txt = txt.replace(
        "} // namespace level23_bh_sync\n",
        "} // namespace level23_bh_sync\n#endif // LEVEL23_BH_SYNC_HELPERS_DEFINED\n",
        1
    )
    age.write_text(txt)
    print("added LEVEL23_BH_SYNC_HELPERS_DEFINED guard to central helper")
PY

cat > run_bigeye_level23_bh_diag_compile_fix_check.sh <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"

echo "== Duplicate helper audit =="
grep -R "namespace level23_bh_sync" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

n=$(grep -R "namespace level23_bh_sync" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' | wc -l | tr -d ' ')

if [[ "$n" != "1" ]]; then
  echo "ERROR: expected exactly one BH helper namespace after cleanup, found $n" >&2
  exit 1
fi

echo
echo "== Required local symbols =="
grep -R "const auto weight = default_weight_at_age\|const auto maturity = default_maturity_at_age\|const double phi0\|expected_recruitment\|next\[0\]" -n \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

echo
echo "== Stale constant recruitment audit =="
if grep -R "next\[0\] = r0\|next\[0\] = r0 \*" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_'; then
  echo "ERROR: stale constant-r0 recruitment remains." >&2
  exit 1
fi

echo
echo "== Build/run Level 23 computed-phi0 BH check =="
if [[ -x ./run_bigeye_level23_unfished_phi0_bh_check.sh ]]; then
  ./run_bigeye_level23_unfished_phi0_bh_check.sh
else
  echo "No run_bigeye_level23_unfished_phi0_bh_check.sh found."
fi
EOS
chmod +x run_bigeye_level23_bh_diag_compile_fix_check.sh

echo
echo "Created: ./run_bigeye_level23_bh_diag_compile_fix_check.sh"
echo "Run:"
echo "  ./run_bigeye_level23_bh_diag_compile_fix_check.sh"
