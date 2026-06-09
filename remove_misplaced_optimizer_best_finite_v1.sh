#!/usr/bin/env bash
set -euo pipefail

echo "== Remove misplaced best-finite line-search return block =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

opt="core/optimizer.hpp"
[[ -f "$opt" ]] || { echo "ERROR: missing $opt" >&2; exit 1; }

cp "$opt" ".quadra_patch_backups/optimizer.hpp.remove_misplaced_best_finite_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("core/optimizer.hpp")
s = p.read_text()

marker = "QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1"
if marker not in s:
    print("Misplaced best-finite block not found; nothing to remove.")
    raise SystemExit(0)

# Remove from comment marker through the matching 'return out;' plus closing brace.
pat = re.compile(
    r'\n\s*// QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1\n'
    r'(?:(?!\n\s*// QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1).)*?'
    r'\n\s*return out;\n\s*}\n',
    re.S,
)

s2, n = pat.subn("\n", s, count=1)
if n != 1:
    raise SystemExit("ERROR: could not remove misplaced block cleanly")

p.write_text(s2)
print("Removed misplaced best-finite block from optimizer.hpp")
PY

cat > inspect_optimizer_misplaced_best_finite_removed_v1.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

echo "== Check misplaced block removed =="
if grep -n "QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1" core/optimizer.hpp; then
  echo "ERROR: misplaced marker still present" >&2
  exit 1
else
  echo "OK: misplaced marker removed"
fi

echo
echo "== Check optimizer still has line-search references =="
grep -n "line search\\|sufficiently decrease\\|runtime_error" core/optimizer.hpp | head -80 || true

echo
echo "Next run:"
echo "  ./inspect_opakapaka_level1_reporting_v7.sh"
SH
chmod +x inspect_optimizer_misplaced_best_finite_removed_v1.sh

echo
echo "Backups saved with suffix: remove_misplaced_best_finite_${stamp}.bak"
echo "Run:"
echo "  ./inspect_optimizer_misplaced_best_finite_removed_v1.sh"
