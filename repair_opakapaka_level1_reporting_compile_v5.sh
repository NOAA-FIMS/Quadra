#!/usr/bin/env bash
set -euo pipefail

echo "== Repair Opakapaka Level-1 reporting compile issues =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
[[ -f "$cpp" ]] || { echo "ERROR: missing $cpp" >&2; exit 1; }

cp "$cpp" ".quadra_patch_backups/opakapaka_projection.cpp.level1_reporting_compile_repair_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/opakapaka_projection/opakapaka_projection.cpp")
s = p.read_text()

# Qualify model types because reporting helpers were inserted outside opakapaka_example.
s = s.replace("std::vector<Observation>& data", "std::vector<opakapaka_example::Observation>& data")
s = s.replace("std::vector<ProjectionRow>& rows", "std::vector<opakapaka_example::ProjectionRow>& rows")

# Repair accidental double-qualification if script is rerun.
s = s.replace("opakapaka_example::opakapaka_example::Observation", "opakapaka_example::Observation")
s = s.replace("opakapaka_example::opakapaka_example::ProjectionRow", "opakapaka_example::ProjectionRow")

# Current source names projected rows 'projection', not 'projection_rows'.
s = s.replace(
    'write_projection_uncertainty_csv("examples/opakapaka_projection/outputs/projection_uncertainty.csv", projection_rows);',
    'write_projection_uncertainty_csv("examples/opakapaka_projection/outputs/projection_uncertainty.csv", projection);'
)

# OptResult contract currently does not expose runtime_ms. Use NaN until runtime is plumbed.
s = s.replace(
    'write_runtime_memory_summary_csv("examples/opakapaka_projection/outputs/runtime_memory_summary.csv", fit.runtime_ms, fit.u_hat.size(), 58);',
    'write_runtime_memory_summary_csv("examples/opakapaka_projection/outputs/runtime_memory_summary.csv", std::numeric_limits<double>::quiet_NaN(), fit.u_hat.size(), 58);'
)

p.write_text(s)
print("Patched compile issues")
PY

cat > inspect_opakapaka_level1_reporting_v5.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
exe="build/examples/opakapaka_projection"

echo "== Relevant reporting lines =="
grep -n "std::vector<opakapaka_example::Observation>\\|std::vector<opakapaka_example::ProjectionRow>\\|projection_uncertainty.csv\\|runtime_memory_summary.csv" "$cpp"

echo
echo "== Build/run Opakapaka example =="

mkdir -p build/examples

eigen_include=""
for d in external/eigen core/eigen external/eigen3 external/Eigen /opt/homebrew/include/eigen3 /usr/local/include/eigen3 /usr/include/eigen3; do
  if [[ -f "$d/Eigen/Core" ]]; then
    eigen_include="$d"
    break
  fi
done

if [[ -z "$eigen_include" ]]; then
  found="$(find . -path '*/Eigen/Core' -type f 2>/dev/null | head -1 || true)"
  [[ -n "$found" ]] && eigen_include="$(dirname "$(dirname "$found")")"
fi

if [[ -z "$eigen_include" ]]; then
  echo "ERROR: could not find Eigen/Core" >&2
  exit 1
fi

echo "Using Eigen include: $eigen_include"

c++ -std=c++17 -O3 -flto \
  -I"$eigen_include" -Icore -I. \
  -o "$exe" "$cpp"

"$exe"

echo
echo "== New outputs =="
ls -1 examples/opakapaka_projection/outputs | grep -E 'uncertainty|covariance|correlation|standard_errors|confidence|derived|runtime' || true

echo
echo "== uncertainty_summary.csv =="
cat examples/opakapaka_projection/outputs/uncertainty_summary.csv

echo
echo "== standard_errors.csv =="
cat examples/opakapaka_projection/outputs/standard_errors.csv

echo
echo "== confidence_intervals.csv =="
cat examples/opakapaka_projection/outputs/confidence_intervals.csv

echo
echo "== derived_quantities.csv preview =="
head -10 examples/opakapaka_projection/outputs/derived_quantities.csv

echo
echo "== projection_uncertainty.csv preview =="
head -10 examples/opakapaka_projection/outputs/projection_uncertainty.csv
SH
chmod +x inspect_opakapaka_level1_reporting_v5.sh

echo
echo "Backups saved with suffix: level1_reporting_compile_repair_${stamp}.bak"
echo "Run:"
echo "  ./inspect_opakapaka_level1_reporting_v5.sh"
