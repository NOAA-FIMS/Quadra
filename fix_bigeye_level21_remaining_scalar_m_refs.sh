#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
TS="$(date +%Y%m%d_%H%M%S)"

FILES=(
  "$ROOT/diagnostics/bigeye_initial_numbers_diagnostics.hpp"
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp"
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"
)

for f in "${FILES[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
  cp "$f" "$f.before_level21_remaining_scalar_m_fix.$TS"
  echo "backup: $f.before_level21_remaining_scalar_m_fix.$TS"
done

python3 - <<'PY'
from pathlib import Path

ROOT = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")

paths = [
    ROOT / "diagnostics/bigeye_initial_numbers_diagnostics.hpp",
    ROOT / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp",
    ROOT / "diagnostics/bigeye_purse_seine_prediction_decomposition.hpp",
    ROOT / "diagnostics/bigeye_longline_prediction_decomposition.hpp",
]

for path in paths:
    text = path.read_text()
    old = text

    text = text.replace(
        'txt << "m_fixed:                    " << m << "\\n";',
        'txt << "m_young:                    " << m_at_age[0] << "\\n";\n'
        '  txt << "m_adult:                    " << m_at_age[3] << "\\n";\n'
        '  txt << "m_old:                      " << m_at_age[kAges - 1] << "\\n";'
    )
    text = text.replace(
        'csv << "summary,m_fixed,," << m << ",\\n";',
        'csv << "summary,m_young,," << m_at_age[0] << ",ages 1-3\\n";\n'
        '  csv << "summary,m_adult,," << m_at_age[3] << ",ages 4-7\\n";\n'
        '  csv << "summary,m_old,," << m_at_age[kAges - 1] << ",ages 8-10\\n";'
    )
    text = text.replace(
        'txt << "m: " << m << "\\n";',
        'txt << "m_young: " << m_at_age[0] << "\\n";\n'
        '  txt << "m_adult: " << m_at_age[3] << "\\n";\n'
        '  txt << "m_old: " << m_at_age[kAges - 1] << "\\n";'
    )

    text = text.replace(
        'const double z_prev = m + fbar * total_sel[prev];',
        'const double z_prev = m_at_age[prev] + fbar * total_sel[prev];'
    )
    text = text.replace(
        'const double z_last = m + fbar * total_sel[last];',
        'const double z_last = m_at_age[last] + fbar * total_sel[last];'
    )

    text = text.replace('        const auto weight = default_weight_at_age();', '    const auto weight = default_weight_at_age();')
    text = text.replace('        const double phi0 =', '    const double phi0 =')
    text = text.replace('\nconst double expected_recruitment =', '\n    const double expected_recruitment =')

    if text != old:
        path.write_text(text)
        print(f"patched: {path}")
    else:
        print(f"NOTE: no changes made to {path}")
PY

cat > run_bigeye_level21_remaining_scalar_m_fix_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Audit undeclared/scalar m remnants in Level 21 diagnostics =="
if grep -R "std::exp(-m)\|1.0 - std::exp(-m)\|z_prev = m +\|z_last = m +\|<< m <<\|m_fixed" -n \
  "$ROOT/diagnostics" \
  --include='*.hpp' --include='*.cpp' \
  | grep -v '\.before_' | grep -v 'before_' >/tmp/level21_remaining_m_refs.txt; then
  echo "WARNING: possible scalar-M remnants remain:"
  cat /tmp/level21_remaining_m_refs.txt
else
  echo "OK: no obvious scalar-M diagnostic remnants."
fi

echo
echo "== Confirm age-specific M dynamics in patched diagnostics =="
grep -R "z_prev = m_at_age\[prev\]\|z_last = m_at_age\[last\]\|m_young:\|summary,m_young" -n \
  "$ROOT/diagnostics" \
  --include='*.hpp' \
  | grep -v '\.before_' | grep -v 'before_' || true

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh
SH

chmod +x run_bigeye_level21_remaining_scalar_m_fix_check.sh

echo
echo "Created:"
echo "  ./run_bigeye_level21_remaining_scalar_m_fix_check.sh"
echo
echo "Run:"
echo "  ./run_bigeye_level21_remaining_scalar_m_fix_check.sh"
