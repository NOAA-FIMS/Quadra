#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
OBJ="$LVL/objective/bigeye_quadra_objective.hpp"
DIAG_DIR="$LVL/diagnostics"
ts="$(date +%Y%m%d_%H%M%S)"

backup() {
  local f="$1"
  [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
  cp "$f" "$f.before_level21_bh_compile_fix_v2.$ts"
  echo "backup: $f.before_level21_bh_compile_fix_v2.$ts"
}

backup "$OBJ"
for f in \
  "$DIAG_DIR/bigeye_age_comp_residual_diagnostics.hpp" \
  "$DIAG_DIR/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$DIAG_DIR/bigeye_longline_prediction_decomposition.hpp" \
  "$DIAG_DIR/bigeye_level21_plus_group_audit.hpp"
do
  [[ -f "$f" ]] && backup "$f"
done

python3 - <<'PY'
from pathlib import Path
import re

lvl = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")
obj = lvl / "objective/bigeye_quadra_objective.hpp"
diag_dir = lvl / "diagnostics"

def read(p): return p.read_text()
def write(p, s): p.write_text(s)

s = read(obj)

if "const auto maturity = default_maturity_at_age();" not in s:
    s = s.replace(
        "const auto weight = default_weight_at_age();",
        "const auto weight = default_weight_at_age();\n    const auto maturity = default_maturity_at_age();",
        1,
    )

if "const T steepness = T(" not in s:
    s = s.replace(
        "const T sigma_rec_dev = T(0.35);",
        "const T sigma_rec_dev = T(0.35);\n    const T steepness = T(0.75);",
        1,
    )

if "const T phi0 =" not in s:
    phi0_block = """    T unfished_spawning_biomass = T(0.0);
    for (int a = 0; a < kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      unfished_spawning_biomass =
          unfished_spawning_biomass + n[i] * T(weight[i]) * T(maturity[i]);
    }
    const T phi0 = unfished_spawning_biomass / max_t(r0, min_positive);

"""
    marker = "    const auto years = unique_years();"
    if marker in s:
        s = s.replace(marker, phi0_block + marker, 1)
        print("inserted objective phi0 before unique_years()")
    else:
        raise SystemExit("ERROR: could not find objective fallback insertion point: const auto years = unique_years();")

write(obj, s)

diag_files = [
    diag_dir / "bigeye_age_comp_residual_diagnostics.hpp",
    diag_dir / "bigeye_purse_seine_prediction_decomposition.hpp",
    diag_dir / "bigeye_longline_prediction_decomposition.hpp",
    diag_dir / "bigeye_level21_plus_group_audit.hpp",
]

for p in diag_files:
    if not p.exists():
        continue
    s = read(p)

    if "level21_bh_sync::beverton_holt_recruitment" not in s and "next[0] = r0 * std::exp(rec_dev);" not in s:
        continue

    s = re.sub(
        r'\n\s+const auto weight = default_weight_at_age\(\);\n\s+const auto maturity = default_maturity_at_age\(\);\n(?=\s*level21_bh_sync::spawning_biomass_from_numbers)',
        "\n",
        s,
    )

    if "const auto weight = default_weight_at_age();" not in s:
        s, n = re.subn(
            r'(const double r0\s*=\s*std::exp\([^;]+;\s*)',
            r'\1\n  const auto weight = default_weight_at_age();\n  const auto maturity = default_maturity_at_age();\n',
            s,
            count=1,
        )
        if n:
            print(f"inserted weight/maturity in {p}")
    elif "const auto maturity = default_maturity_at_age();" not in s:
        s = s.replace(
            "const auto weight = default_weight_at_age();",
            "const auto weight = default_weight_at_age();\n  const auto maturity = default_maturity_at_age();",
            1,
        )
        print(f"inserted maturity in {p}")

    if "const double phi0 =" not in s:
        s2, n = re.subn(
            r'(std::array<double,\s*kAges>\s+n\s*=\s*unfished_equilibrium_numbers\([^;]+;\s*)',
            r'\1\n  const double phi0 =\n      level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity) /\n      std::max(r0, level21_bh_sync::bh_min_positive());\n',
            s,
            count=1,
            flags=re.S,
        )
        if n:
            s = s2
            print(f"inserted phi0 in {p}")

    if "next[0] = r0 * std::exp(rec_dev);" in s:
        s = s.replace(
            "next[0] = r0 * std::exp(rec_dev);",
            """const double spawning_biomass =
        level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);
    const double expected_recruitment =
        level21_bh_sync::beverton_holt_recruitment(spawning_biomass, r0, 0.75, phi0);
    next[0] = expected_recruitment * std::exp(rec_dev);""",
            1,
        )
        print(f"patched stale recruitment in {p}")

    write(p, s)

print("done")
PY

cat > run_bigeye_level21_bh_compile_fix_v2_check.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Objective BH/phi0 context =="
grep -n "maturity\|steepness\|unfished_spawning_biomass\|phi0\|spawning_biomass\|expected_recruitment\|next\[0\]" \
  "$LVL/objective/bigeye_quadra_objective.hpp" | head -140

echo
echo "== Stale constant recruitment audit =="
if grep -R "next\[0\] = r0\|next\[0\] = r0 \* std::exp\|next\[0\] = r0 \* exp_t" -n \
  "$LVL" --include='*.hpp' --include='*.cpp' | grep -v '\.before_'; then
  echo "ERROR: stale constant recruitment remains."
  exit 1
else
  echo "OK: no stale constant-R0 recruitment references found."
fi

echo
echo "== Build/run Level 21 BH check =="
./run_bigeye_level21_age_based_m_check.sh
EOF

chmod +x run_bigeye_level21_bh_compile_fix_v2_check.sh

echo
echo "Created:"
echo "  ./run_bigeye_level21_bh_compile_fix_v2_check.sh"
echo
echo "Run:"
echo "  ./run_bigeye_level21_bh_compile_fix_v2_check.sh"
