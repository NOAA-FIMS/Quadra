#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/afsc_walleye_pollock"
QUADRA="$BASE/quadra"
CPP="$QUADRA/walleye_pollock.cpp"
MODEL="$QUADRA/model/pollock_model.hpp"

if [[ ! -f "$CPP" ]]; then
  echo "ERROR: $CPP not found. Run from repo root."
  exit 1
fi

mkdir -p "$QUADRA/model"

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$CPP" "${CPP}.before_recover_pollock_model_py39.${STAMP}"
[[ -f "$MODEL" ]] && cp "$MODEL" "${MODEL}.before_recover_pollock_model_py39.${STAMP}" || true

python3 - <<'PYCODE'
from pathlib import Path
import sys

base = Path("examples/NMFS/afsc_walleye_pollock")
quadra = base / "quadra"
cpp = quadra / "walleye_pollock.cpp"
model = quadra / "model" / "pollock_model.hpp"

def find_block(text, signature):
    start = text.find(signature)
    if start == -1:
        return None
    brace = text.find("{", start)
    if brace == -1:
        return None
    depth = 0
    end = None
    for i in range(brace, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                semi = text.find(";", i)
                if semi == -1:
                    return None
                end = semi + 1
                break
    if end is None:
        return None
    block = text[start:end]
    lines = []
    for line in block.splitlines():
        lines.append(line[2:] if line.startswith("  ") else line)
    return "\n".join(lines)

candidates = [cpp]
candidates += sorted(quadra.glob("walleye_pollock.cpp.before*"), reverse=True)
candidates += sorted(Path("patch_backups").glob("**/walleye_pollock.cpp.before*"), reverse=True)

obs_block = None
model_block = None
source = None

for cand in candidates:
    try:
        text = cand.read_text()
    except Exception:
        continue

    obs = find_block(text, "  struct Obs")
    if obs is None:
        obs = find_block(text, "struct Obs")

    pm = find_block(text, "  struct PollockModel")
    if pm is None:
        pm = find_block(text, "struct PollockModel")

    if obs and pm:
        obs_block = obs
        model_block = pm
        source = cand
        break

if not obs_block or not model_block:
    print("ERROR: could not find a backup containing both struct Obs and struct PollockModel.")
    print("Searched first candidates:")
    for c in candidates[:30]:
        print(" ", c)
    sys.exit(1)

old_constants = """      constexpr int A = 7;
      const double weight[A] = {0.20, 0.45, 0.75, 1.10, 1.45, 1.75, 2.00};
      const double maturity[A] = {0.00, 0.10, 0.45, 0.80, 0.95, 1.00, 1.00};
      const double M = 0.25;

"""
model_block = model_block.replace(old_constants, "")

repls = {
    "std::vector<AD> N(A);": "std::vector<AD> N(pollock::n_ages);",
    "std::vector<AD> caa(A);": "std::vector<AD> caa(pollock::n_ages);",
    "std::vector<AD> next(A);": "std::vector<AD> next(pollock::n_ages);",
    "N[A - 2]": "N[pollock::n_ages - 2]",
    "A - 1": "pollock::n_ages - 1",
    "a < A": "a < pollock::n_ages",
    "AD(M)": "AD(pollock::natural_mortality)",
    "weight[a]": "pollock::weight_at_age[a]",
    "maturity[a]": "pollock::maturity_at_age[a]",
}
for old, new in repls.items():
    model_block = model_block.replace(old, new)

header = """#pragma once

#include "pollock_constants.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

"""

model.write_text(header + obs_block + "\n\n" + model_block + "\n")

s = cpp.read_text()

for stale in [
    '#include "pollock_model.hpp"\n',
    '#include "../model/pollock_model.hpp"\n',
]:
    s = s.replace(stale, "")

if '#include "model/pollock_model.hpp"\n' not in s:
    if '#include "model/pollock_constants.hpp"\n' in s:
        s = s.replace(
            '#include "model/pollock_constants.hpp"\n',
            '#include "model/pollock_constants.hpp"\n#include "model/pollock_model.hpp"\n',
            1,
        )
    else:
        lines = s.splitlines()
        insert_at = 0
        while insert_at < len(lines) and lines[insert_at].startswith("#include"):
            insert_at += 1
        lines.insert(insert_at, '#include "model/pollock_model.hpp"')
        s = "\n".join(lines) + "\n"

s = s.replace("using pollock_example::Obs;\n", "")
s = s.replace("using pollock_example::PollockModel;\n", "")

def remove_struct_block(text, signature):
    while True:
        start = text.find(signature)
        if start == -1:
            return text
        brace = text.find("{", start)
        semi_before_brace = text.find(";", start, brace if brace != -1 else start)
        if brace == -1 or (semi_before_brace != -1 and semi_before_brace < brace):
            return text
        depth = 0
        end = None
        for i in range(brace, len(text)):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
                if depth == 0:
                    semi = text.find(";", i)
                    if semi == -1:
                        return text
                    end = semi + 1
                    break
        if end is None:
            return text
        while end < len(text) and text[end] in " \t\r\n":
            end += 1
        text = text[:start] + "\n\n" + text[end:]

for sig in ["  struct Obs", "struct Obs", "  struct PollockModel", "struct PollockModel"]:
    s = remove_struct_block(s, sig)

lines = s.splitlines()
seen = set()
out = []
for line in lines:
    if line.startswith("#include "):
        if line in seen:
            continue
        seen.add(line)
    out.append(line)
cpp.write_text("\n".join(out) + "\n")

print("Recovered Obs and PollockModel from: {}".format(source))
PYCODE

cat > inspect_recovered_pollock_model.sh <<'EOF_INSPECT'
#!/usr/bin/env bash
set -euo pipefail

echo "== Active cpp top =="
sed -n '1,28p' examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp

echo
echo "== Model header declarations =="
grep -n "struct Obs\|struct PollockModel\|namespace pollock_example" \
  examples/NMFS/afsc_walleye_pollock/quadra/model/pollock_model.hpp || true

echo
echo "== Model header first 80 lines =="
sed -n '1,80p' examples/NMFS/afsc_walleye_pollock/quadra/model/pollock_model.hpp

echo
echo "== Local type definition checks in cpp =="
if grep -n "struct Obs\|struct PollockModel" \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp; then
  echo "WARNING: active cpp still has local type definitions"
else
  echo "OK: active cpp has no local Obs/PollockModel definitions"
fi
EOF_INSPECT

chmod +x inspect_recovered_pollock_model.sh

cat > run_recovered_pollock_model_check.sh <<'EOF_RUN'
#!/usr/bin/env bash
set -euo pipefail

./inspect_recovered_pollock_model.sh

echo
echo "== Rebuild O3 Pollock showcase after recovering model header =="
./run_pollock_driver_showcase_report.sh
EOF_RUN

chmod +x run_recovered_pollock_model_check.sh

echo "Recovered real Obs/PollockModel into quadra/model/pollock_model.hpp."
echo
echo "Run:"
echo "  ./inspect_recovered_pollock_model.sh"
echo "  ./run_recovered_pollock_model_check.sh"
