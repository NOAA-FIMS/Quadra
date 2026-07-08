#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp"
DIR="$(dirname "$OBJ")"

echo "== Current first 30 lines =="
sed -n '1,30p' "$OBJ" || true

echo
echo "== Candidate backups =="
ls -1t "$DIR"/bigeye_quadra_objective.hpp* 2>/dev/null || true

echo
echo "== Inspect candidate headers =="
python3 - <<'PY'
from pathlib import Path

obj = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp")
candidates = sorted(obj.parent.glob(obj.name + "*"), key=lambda p: p.stat().st_mtime, reverse=True)

good = []
for p in candidates:
    text = p.read_text(errors="replace")
    first = "\n".join(text.splitlines()[:8])
    is_good = (
        "#pragma once" in text[:500]
        and "class BigeyeQuadraObjective" in text
        and "template <class T>" in text
        and "T operator()" in text
    )
    print(f"\n--- {p} ---")
    print(f"good_candidate={is_good}")
    print(first)
    if is_good and p != obj:
        good.append(p)

if not good:
    raise SystemExit("\nERROR: no good backup found. Need restore from git or copy Level 21 objective from a prior saved file.")

restore = good[0]
obj.write_text(restore.read_text())
print(f"\nRESTORED_FROM={restore}")
PY

echo
echo "== Restored first 30 lines =="
sed -n '1,30p' "$OBJ"

echo
echo "== Compile smoke test after restore =="
./run_bigeye_level21_age_based_m_check.sh
