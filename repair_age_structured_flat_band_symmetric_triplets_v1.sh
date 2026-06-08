#!/usr/bin/env bash
set -euo pipefail

target="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_age_structured_no_plus_flat_band_benchmark_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/benchmark_age_structured_no_plus_flat_band.cpp.symmetric.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp")
s = p.read_text()

old = '''      const double v = Hband[static_cast<std::size_t>((d + bw) * n + i)];
      if (std::abs(v) > 1e-12) {
        triplets.emplace_back(i, j, v);
      }
'''

new = '''      const double v = Hband[static_cast<std::size_t>((d + bw) * n + i)];
      if (std::abs(v) > 1e-12) {
        triplets.emplace_back(i, j, v);
        if (i != j) {
          triplets.emplace_back(j, i, v);
        }
      }
'''

if old not in s:
    raise SystemExit("Could not find flat-band triplet emission block")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Fixed flat-band symmetric triplet emission.

Run:
  ./run_quadra_flat_band_vs_tmb_age_structured_no_plus_benchmark.sh 10 25,50,100,250,500,1000 10

Expected:
  objectives match TMB
  nnz roughly doubles from the bad flat-band result

EOF
