#!/usr/bin/env bash
set -euo pipefail

file="core/laplace.hpp"
backup="${file}.backup_force_remove_bad_laplace_tail_block_$(date +%Y%m%d_%H%M%S)"

if [[ ! -f "$file" ]]; then
  echo "ERROR: $file not found. Run from Quadra repo root." >&2
  exit 1
fi

cp "$file" "$backup"

python3 - <<'PY'
from pathlib import Path
p = Path('core/laplace.hpp')
s = p.read_text()

# Remove the accidentally inserted tail diagnostic block.  The compile error
# shows it contains timing_hdot_end/timing_hdot_start and ends with return grad;
# in a scope where grad does not exist.
idx = s.find('timing_hdot_end - timing_hdot_start')
if idx == -1:
    print('No timing_hdot_end block found; file may already be clean.')
    raise SystemExit(0)

# Walk backward to the start of the diagnostic/logging statement/block.
# Prefer a nearby preprocessor guard if present, otherwise the nearest cerr line.
window_start = max(0, idx - 2500)
prefix = s[window_start:idx]
starts = []
for marker in ['#ifdef QUADRA_GRADIENT_DIAGNOSTIC', '#if defined(QUADRA_GRADIENT_DIAGNOSTIC)', 'std::cerr']:
    j = prefix.rfind(marker)
    if j != -1:
        starts.append(window_start + j)
if not starts:
    raise RuntimeError('Found timing_hdot_end but could not identify block start.')
start = min(starts) if any(s[t:t+1] == '#' for t in starts) else max(starts)

# Walk forward through the bad return grad; line.  This removes the whole bad tail.
ret = s.find('return grad;', idx)
if ret == -1:
    raise RuntimeError('Found timing_hdot_end but not following return grad;')
line_end = s.find('\n', ret)
if line_end == -1:
    line_end = len(s)
else:
    line_end += 1

removed = s[start:line_end]
new = s[:start] + s[line_end:]
p.write_text(new)

print('Removed bad block from core/laplace.hpp')
print('Removed lines containing:')
for line in removed.splitlines():
    if 'timing_hdot' in line or 'return grad' in line or 'QUADRA_GRADIENT_DIAGNOSTIC' in line:
        print('  ' + line.strip())
PY

echo "Backup saved to: $backup"
echo "Now run:"
echo "  grep -n \"timing_hdot_end\\|return grad;\" core/laplace.hpp"
echo "Then rebuild."
