#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
OBJ="$LVL/objective/bigeye_quadra_objective.hpp"
DIAG_DIR="$LVL/diagnostics"
ts="$(date +%Y%m%d_%H%M%S)"

backup() {
  local f="$1"
  [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
  cp "$f" "$f.before_level21_bh_compile_fix.$ts"
  echo "backup: $f.before_level21_bh_compile_fix.$ts"
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
def write(p,s): p.write_text(s)

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
    patterns = [
        r'(std::array<T,\s*kAges>\s+n\s*=\s*unfished_equilibrium_numbers\([^;]+;\s*)',
        r'(auto\s+n\s*=\s*unfished_equilibrium_numbers\([^;]+;\s*)',
    ]
    inserted = False
    block = r'''\1
    T unfished_spawning_biomass = T(0.0);
    for (int a = 0; a < kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      unfished_spawning_biomass =
          unfished_spawning_biomass + n[i] * T(weight[i]) * T(maturity[i]);
    }
    const T phi0 = unfished_spawning_biomass / max_t(r0, min_positive);
'''
    for pat in patterns:
        s2, n = re.subn(pat, block, s, count=1, flags=re.S)
        if n == 1:
            s = s2
            inserted = True
            print("inserted objective phi0")
            break
    if not inserted:
        raise SystemExit("ERROR: could not find objective insertion point for phi0")

write(obj, s)

for p in [
    diag_dir / "bigeye_age_comp_residual_diagnostics.hpp",
    diag_dir / "bigeye_purse_seine_prediction_decomposition.hpp",
    diag_dir / "bigeye_longline_prediction_decomposition.hpp",
]:
    if not p.exists():
        continue
    s = read(p)
    if "level21_bh_sync::beverton_holt_recruitment" not in s:
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
        if n == 0:
            print(f"WARNING: no r0 insertion point for weight/maturity in {p}")

    if "const auto weight = default_weight_at_age();" in s and "const auto maturity = default_maturity_at_age();" not in s:
        s = s.replace(
            "const auto weight = default_weight_at_age();",
            "const auto weight = default_weight_at_age();\n  const auto maturity = default_maturity_at_age();",
            1,
        )

    if "const double phi0 =" not in s:
        s, n = re.subn(
            r'(std::array<double,\s*kAges>\s+n\s*=\s*unfished_equilibrium_numbers\([^;]+;\s*)',
            r'\1\n  const double phi0 =\n      level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity) /\n      std::max(r0, level21_bh_sync::bh_min_positive());\n',
            s,
            count=1,
            flags=re.S,
        )
        if n == 0:
            print(f"WARNING: no n insertion point for phi0 in {p}")

    write(p, s)

p = diag_dir / "bigeye_level21_plus_group_audit.hpp"
if p.exists():
    s = read(p)
    if "next[0] = r0 * std::exp(rec_dev);" in s:
        if "const auto maturity = default_maturity_at_age();" not in s and "const auto weight = default_weight_at_age();" in s:
            s = s.replace(
                "const auto weight = default_weight_at_age();",
                "const auto weight = default_weight_at_age();\n  const auto maturity = default_maturity_at_age();",
                1,
            )

        if "const double phi0 =" not in s:
            s, n = re.subn(
                r'(std::array<double,\s*kAges>\s+n\s*=\s*unfished_equilibrium_numbers\([^;]+;\s*)',
                r'\1\n  const double phi0 =\n      level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity) /\n      std::max(r0, level21_bh_sync::bh_min_positive());\n',
                s,
                count=1,
                flags=re.S,
            )
            if n == 0:
                print("WARNING: no plus-group audit n insertion point for phi0")

        s = s.replace(
            "next[0] = r0 * std::exp(rec_dev);",
            '''const double spawning_biomass =
        level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);
    const double expected_recruitment =
        level21_bh_sync::beverton_holt_recruitment(spawning_biomass, r0, 0.75, phi0);
    next[0] = expected_recruitment * std::exp(rec_dev);''',
            1,
        )
        print("patched stale plus-group audit recruitment")
    write(p, s)

print("done")
PY

cat > run_bigeye_level21_bh_compile_fix_check.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Objective BH/phi0 context =="
grep -n "maturity\|steepness\|unfished_spawning_biomass\|phi0\|spawning_biomass\|expected_recruitment\|next\[0\]" \
  "$LVL/objective/bigeye_quadra_objective.hpp" | head -120

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

chmod +x run_bigeye_level21_bh_compile_fix_check.sh

echo
echo "Created:"
echo "  ./run_bigeye_level21_bh_compile_fix_check.sh"
echo
echo "Run:"
echo "  ./run_bigeye_level21_bh_compile_fix_check.sh"
