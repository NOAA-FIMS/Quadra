#!/usr/bin/env bash
set -euo pipefail

# repair_structure_detector_banded_test_threshold_v1.sh
#
# Fixes the banded registry test.
#
# The test constructed n=20, bandwidth=4, which creates about 160 nonzeros:
#   fill ~= 160 / 400 = 0.40+
#
# With dense_fill_ratio = 0.40, the detector correctly avoids the banded
# backend and falls toward dense/sparse behavior. The test expectation was
# too strict for its own threshold.
#
# This patch increases the test size to n=40 for the same bandwidth, dropping
# fill ratio well below the dense threshold while preserving the banded pattern.

target="tests/test_structure_detector_registry.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_structure_detector_registry_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/test_structure_detector_registry.cpp.banded_threshold.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

old = """void test_banded_recommendation() {
  const int n = 20;
  std::vector<Eigen::Triplet<double>> t;"""

new = """void test_banded_recommendation() {
  const int n = 40;
  std::vector<Eigen::Triplet<double>> t;"""

if old not in s:
    raise SystemExit("Could not find banded test n=20 block")

s = s.replace(old, new, 1)
p.write_text(s)
PY

cat <<'EOF'

Repaired banded detector test threshold.

Run:
  ./run_structure_detector_registry_test.sh

EOF
