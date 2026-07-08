#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna"
STAMP="$(date +%Y%m%d_%H%M%S)"

targets=(
  "$ROOT/level16_purse_seine_age_selectivity/quadra/bigeye_age_structured.hpp"
  "$ROOT/level17_juvenile_mortality_diagnostic/quadra/bigeye_age_structured.hpp"
)

for f in "${targets[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "missing: $f"
    continue
  fi

  cp "$f" "${f}.before_age_comp_reader_fix.${STAMP}"

  python3 - "$f" <<'PY'
from pathlib import Path
import sys

p = Path(sys.argv[1])
s = p.read_text()

# CSV header is:
# year,fleet,catch_mt,index,age_comp_1,...,age_comp_10
# That is 14 fields, with age comps starting at field index 4.
s = s.replace(
    'if (fields.size() != 13) {\n      throw std::runtime_error("Expected 13 columns in observations CSV");\n    }',
    'if (fields.size() != 14) {\n      throw std::runtime_error("Expected 14 columns in observations CSV");\n    }'
)

# Old single-fleet reader expected year,catch,index,age comps...
# For fleet CSVs, field[1] is fleet, field[2] is catch, field[3] is index,
# and age comps start at field[4].
s = s.replace(
    'obs.catch_mt = std::stod(fields[1]);\n    obs.index = std::stod(fields[2]);',
    'obs.catch_mt = std::stod(fields[2]);\n    obs.index = std::stod(fields[3]);'
)

s = s.replace(
    'obs.age_comp[static_cast<std::size_t>(a)] = std::stod(fields[3 + a]);',
    'obs.age_comp[static_cast<std::size_t>(a)] = std::stod(fields[4 + a]);'
)

# If the struct has no fleet field, leave it alone; this helper may only be
# used for generic single-observation reports. The important part is not
# polluting age_comp_1 with the index value.

p.write_text(s)
PY

  echo "patched: $f"
done

echo
echo "== remaining stale field readers =="
grep -R "fields\\[3 + a\\]\\|Expected 13 columns" -n \
  "$ROOT/level16_purse_seine_age_selectivity" \
  "$ROOT/level17_juvenile_mortality_diagnostic" || true

echo
echo "Now rerun Level 16 or Level 17 checks if you want to confirm no behavior changed in the active driver path."
