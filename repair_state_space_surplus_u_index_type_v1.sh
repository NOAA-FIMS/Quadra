#!/usr/bin/env bash
set -euo pipefail

# repair_state_space_surplus_u_index_type_v1.sh
#
# Fixes compile error:
#   no matching function for call to 'max'
# caused by mixing int and Eigen::Index in:
#   std::max(1, u.size() - 1)

target="examples/state_space_surplus_production/fit_state_space_surplus_u.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_state_space_surplus_optimize_u_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/fit_state_space_surplus_u.cpp.index_type.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("examples/state_space_surplus_production/fit_state_space_surplus_u.cpp")
s = p.read_text()

old = "  sd = std::sqrt(sd / static_cast<double>(std::max(1, u.size() - 1)));"
new = """  const int sd_denominator = std::max(1, static_cast<int>(u.size()) - 1);
  sd = std::sqrt(sd / static_cast<double>(sd_denominator));"""

if old not in s:
    raise SystemExit("Could not find std::max denominator line to patch")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Fixed Eigen::Index/std::max type mismatch.

Run:
  ./run_fit_state_space_surplus_u_example.sh

EOF
