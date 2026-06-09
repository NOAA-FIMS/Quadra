#!/usr/bin/env bash
set -euo pipefail

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
exe="build/examples/opakapaka_projection"
ad_impl="build/examples/opakapaka_adgraph_global.cpp"

echo "== H_uu current API repair markers =="
grep -n "HUU_CURRENT_API_REPAIR\\|ADScope scope(graph)\\|p_full.emplace_back\\|scope.backward\\|get_pattern\\|extract_sparse_hessian(scope" "$cpp" || true

bad=0
grep -q "scope.variable" "$cpp" && { echo "ERROR: stale scope.variable remains" >&2; bad=1; }
grep -q "discover_sparse_hessian_pattern" "$cpp" && { echo "ERROR: stale discover_sparse_hessian_pattern remains" >&2; bad=1; }
grep -q "PropagateAdjoint(objective)" "$cpp" && { echo "ERROR: stale PropagateAdjoint(objective) remains" >&2; bad=1; }
[[ "$bad" -eq 0 ]] || exit 1

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
echo "== random_effect_uncertainty.csv preview =="
head -12 examples/opakapaka_projection/outputs/random_effect_uncertainty.csv

echo
echo "== Check conditional SEs populated =="
if grep -q "pending selected-inverse" examples/opakapaka_projection/outputs/random_effect_uncertainty.csv; then
  echo "ERROR: pending selected-inverse text still present" >&2
  exit 1
fi

if awk -F, 'NR > 1 && ($3 == "" || $3 == "nan" || $3 == "NaN") { bad=1 } END { exit bad }' \
  examples/opakapaka_projection/outputs/random_effect_uncertainty.csv; then
  echo "OK: random-effect conditional SEs are populated"
else
  echo "ERROR: at least one random-effect conditional SE is missing/NaN" >&2
  exit 1
fi
