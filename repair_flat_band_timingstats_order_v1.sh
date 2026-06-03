#!/usr/bin/env bash
set -euo pipefail

target="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_timing.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_age_structured_flat_band_timing_diagnostics_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/benchmark_age_structured_no_plus_flat_band_timing.cpp.timingstats.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

timing_struct = """struct TimingStats {
  int newton_iterations = 0;
  int eval_all_calls = 0;
  int line_search_eval_calls = 0;
};

"""

# Remove existing TimingStats block wherever it landed.
s = s.replace(timing_struct, "")

marker = "Eigen::VectorXd optimize_x(const Data& data, const Parameters& par, TimingStats* stats = nullptr) {"
if marker not in s:
    raise SystemExit("Could not find optimize_x with TimingStats")

s = s.replace(marker, timing_struct + marker, 1)

p.write_text(s)
PY

cat <<'EOF'

Fixed TimingStats declaration order.

Run:
  ./run_quadra_age_structured_no_plus_flat_band_timing.sh 5 25,50,100,250,500,1000 10

EOF
