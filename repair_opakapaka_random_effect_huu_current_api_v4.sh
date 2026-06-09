#!/usr/bin/env bash
set -euo pipefail

echo "== Repair Opakapaka H_uu extraction to match current ADScope API =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
[[ -f "$cpp" ]] || { echo "ERROR: missing $cpp" >&2; exit 1; }

cp "$cpp" ".quadra_patch_backups/opakapaka_projection.cpp.huu_current_api_repair_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/opakapaka_projection/opakapaka_projection.cpp")
s = p.read_text()

if "QUADRA_OPAKAPAKA_HUU_CURRENT_API_REPAIR_V1" in s:
    print("Current ADScope API repair already installed.")
    raise SystemExit(0)

old = r'''  quadra::ADScope scope;
  std::vector<quadra::AD> p_full;
  p_full.reserve(n_total);

  for (std::size_t i = 0; i < n_fixed; ++i) {
    p_full.push_back(scope.variable(fit.par.at(i)));
  }
  for (std::size_t i = 0; i < n_random; ++i) {
    p_full.push_back(scope.variable(fit.u_hat.at(i)));
  }

  had::g_ADGraph = &scope.graph;

  const auto objective = model(p_full);
  scope.graph.Forward();
  had::PropagateAdjoint(objective);

  auto pattern = quadra::discover_sparse_hessian_pattern(scope, p_full, random_idx);
  auto h_uu = quadra::extract_sparse_hessian(scope, p_full, random_idx, pattern);
'''

new = r'''  // QUADRA_OPAKAPAKA_HUU_CURRENT_API_REPAIR_V1
  had::ADGraph graph;
  quadra::ADScope scope(graph);

  std::vector<quadra::AD> p_full;
  p_full.reserve(n_total);

  for (std::size_t i = 0; i < n_fixed; ++i) {
    p_full.emplace_back(quadra::AD(fit.par.at(i)));
  }
  for (std::size_t i = 0; i < n_random; ++i) {
    p_full.emplace_back(quadra::AD(fit.u_hat.at(i)));
  }

  quadra::AD nll = model(p_full);
  scope.backward(nll);

  const auto& pattern = quadra::get_pattern(scope, p_full, random_idx);
  auto h_uu = quadra::extract_sparse_hessian(scope, p_full, random_idx, pattern);
'''

if old in s:
    s = s.replace(old, new, 1)
else:
    s = re.sub(r'quadra::ADScope\s+scope\s*;', 'had::ADGraph graph;\n  quadra::ADScope scope(graph);', s, count=1)
    s = s.replace('p_full.push_back(scope.variable(fit.par.at(i)));',
                  'p_full.emplace_back(quadra::AD(fit.par.at(i)));')
    s = s.replace('p_full.push_back(scope.variable(fit.u_hat.at(i)));',
                  'p_full.emplace_back(quadra::AD(fit.u_hat.at(i)));')
    s = s.replace('const auto objective = model(p_full);\n  scope.graph.Forward();\n  had::PropagateAdjoint(objective);',
                  'quadra::AD nll = model(p_full);\n  scope.backward(nll);')
    s = s.replace('auto pattern = quadra::discover_sparse_hessian_pattern(scope, p_full, random_idx);',
                  'const auto& pattern = quadra::get_pattern(scope, p_full, random_idx);')
    s = s.replace('had::ADGraph graph;\n  quadra::ADScope scope(graph);',
                  '// QUADRA_OPAKAPAKA_HUU_CURRENT_API_REPAIR_V1\n  had::ADGraph graph;\n  quadra::ADScope scope(graph);',
                  1)

p.write_text(s)
print("Patched H_uu extraction to current ADScope API")
PY

cat > inspect_opakapaka_random_effect_selected_inverse_v4.sh <<'SH'
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
SH
chmod +x inspect_opakapaka_random_effect_selected_inverse_v4.sh

echo
echo "Backups saved with suffix: huu_current_api_repair_${stamp}.bak"
echo "Run:"
echo "  ./inspect_opakapaka_random_effect_selected_inverse_v4.sh"
