#!/usr/bin/env bash
set -euo pipefail

FILE="core/optimizer.hpp"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run from repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${FILE}.before_fix_lbfgs_recovery_compile.${STAMP}"
cp "$FILE" "$BACKUP"
echo "Backed up $FILE to:"
echo "  $BACKUP"

python3 - <<'PY'
from pathlib import Path

p = Path("core/optimizer.hpp")
lines = p.read_text().splitlines()
out = []
i = 0

while i < len(lines):
    line = lines[i]

    # Repair split stream string literals from the previous patch.
    if line.rstrip() == '        std::cout << "L-BFGS runtime_error caught: " << msg << "':
        out.append('        std::cout << "L-BFGS runtime_error caught: " << msg << "\\n";')
        if i + 1 < len(lines) and lines[i + 1].strip() == '";':
            i += 2
            continue

    if line.rstrip() == '        std::cout << "  gnorm = " << gnorm << "':
        out.append('        std::cout << "  gnorm = " << gnorm << "\\n";')
        if i + 1 < len(lines) and lines[i + 1].strip() == '";':
            i += 2
            continue

    if line.rstrip() == '        std::cout << "  max|grad| = " << max_grad << "':
        out.append('        std::cout << "  max|grad| = " << max_grad << "\\n";')
        if i + 1 < len(lines) and lines[i + 1].strip() == '";':
            i += 2
            continue

    # Remove non-assignable LBFGSSolver restart assignment.
    if line.strip() == 'LBFGSSolver<double> restarted_solver(param);':
        if i + 1 < len(lines) and lines[i + 1].strip() == 'solver = std::move(restarted_solver);':
            out.append('            // LBFGSSolver stores param by reference and is not assignable.')
            out.append('            // Calling minimize() again from the accepted recovery point')
            out.append('            // rebuilds the quasi-Newton history inside LBFGSpp.')
            i += 2
            continue

    out.append(line)
    i += 1

p.write_text("\n".join(out) + "\n")
print("Fixed L-BFGS recovery compile issues.")
PY

echo
echo "Rebuild:"
echo 'clang++ -std=c++17 -g -I"external/eigen/" \'
echo '  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_projection.cpp \'
echo '  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \'
echo '  -o build/examples/pifsc_opakapaka'
