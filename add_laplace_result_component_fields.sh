#!/usr/bin/env bash
set -euo pipefail

FILE="core/laplace.hpp"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run this from the Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${FILE}.before_laplace_result_component_fields.${STAMP}"
cp "$FILE" "$BACKUP"
echo "Backed up $FILE to:"
echo "  $BACKUP"

python3 - <<'PY'
from pathlib import Path
import re

path = Path("core/laplace.hpp")
text = path.read_text()

if "joint_objective" in text and "laplace_logdet" in text and "laplace_constant" in text:
    print("LaplaceResult component fields already appear to exist. No patch needed.")
    raise SystemExit(0)

m = re.search(
    r'(template\s*<\s*typename\s+Model\s*>\s*\n\s*struct\s+LaplaceResult\s*\{)',
    text
)
if not m:
    m = re.search(r'(struct\s+LaplaceResult\s*\{)', text)

if not m:
    raise RuntimeError("Could not find struct LaplaceResult in core/laplace.hpp")

insert_at = m.end()

fields = """\n
  // Component breakdown of the Laplace objective:
  //
  //   value = joint_objective + 0.5 * laplace_logdet - laplace_constant
  //
  // These are intentionally stored for diagnostics/reporting and for
  // optimizer-side bookkeeping. They do not change the objective math.
  double joint_objective = 0.0;
  double laplace_logdet = 0.0;
  double laplace_constant = 0.0;
"""

text = text[:insert_at] + fields + text[insert_at:]
path.write_text(text)
print("Inserted joint_objective, laplace_logdet, and laplace_constant into LaplaceResult.")
PY

echo
echo "Relevant LaplaceResult region:"
grep -n "struct LaplaceResult\|joint_objective\|laplace_logdet\|laplace_constant" "$FILE" | head -40

echo
echo "Done. Rebuild now:"
echo 'clang++ -std=c++17 -g -I"external/eigen/" examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp examples/NMFS/sefsc_red_snapper/quadra/red_snapper_age_structured.cpp'
