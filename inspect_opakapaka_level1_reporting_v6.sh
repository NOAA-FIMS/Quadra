#!/usr/bin/env bash
set -euo pipefail

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
exe="build/examples/opakapaka_projection"
ad_impl="build/examples/opakapaka_adgraph_global.cpp"

echo "== Relevant reporting lines =="
grep -n "QUADRA_LEVEL1_UNCERTAINTY_REPORTING_V3\\|write_uncertainty_summary_csv\\|write_projection_uncertainty_csv" "$cpp"

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

cat > "$ad_impl" <<'CPP'
// Standalone diagnostic/example build support.
// Quadra's had_quadra.hpp declares had::g_ADGraph as a thread-local pointer.
#include "core/had_quadra.hpp"

namespace had {
threadDefine ADGraph* g_ADGraph = nullptr;
}
CPP

c++ -std=c++17 -O3 -flto \
  -I"$eigen_include" -Icore -I. \
  -o "$exe" "$cpp" "$ad_impl"

"$exe"

echo
echo "== New outputs =="
ls -1 examples/opakapaka_projection/outputs | grep -E 'uncertainty|covariance|correlation|standard_errors|confidence|derived|runtime' || true

echo
echo "== uncertainty_summary.csv =="
cat examples/opakapaka_projection/outputs/uncertainty_summary.csv

echo
echo "== covariance_matrix.csv =="
cat examples/opakapaka_projection/outputs/covariance_matrix.csv

echo
echo "== standard_errors.csv =="
cat examples/opakapaka_projection/outputs/standard_errors.csv

echo
echo "== confidence_intervals.csv =="
cat examples/opakapaka_projection/outputs/confidence_intervals.csv

echo
echo "== random_effect_uncertainty.csv preview =="
head -10 examples/opakapaka_projection/outputs/random_effect_uncertainty.csv

echo
echo "== derived_quantities.csv preview =="
head -10 examples/opakapaka_projection/outputs/derived_quantities.csv

echo
echo "== projection_uncertainty.csv preview =="
head -10 examples/opakapaka_projection/outputs/projection_uncertainty.csv
