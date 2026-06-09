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
