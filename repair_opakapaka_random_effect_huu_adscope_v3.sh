#!/usr/bin/env bash
set -euo pipefail

echo "== Repair Opakapaka H_uu extraction using ADScope + AD vector API =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
lap="core/laplace.hpp"
[[ -f "$cpp" ]] || { echo "ERROR: missing $cpp" >&2; exit 1; }
[[ -f "$lap" ]] || { echo "ERROR: missing $lap" >&2; exit 1; }

cp "$cpp" ".quadra_patch_backups/opakapaka_projection.cpp.huu_adscope_repair_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/opakapaka_projection/opakapaka_projection.cpp")
s = p.read_text()

if "QUADRA_OPAKAPAKA_HUU_ADSCOPE_REPAIR_V1" in s:
    print("ADScope H_uu repair already installed.")
    raise SystemExit(0)

# Replace entire compute_final_random_effect_hessian function body/signature block.
pat = re.compile(
    r'template\s*<class\s+Model>\s*'
    r'Eigen::SparseMatrix<double>\s+compute_final_random_effect_hessian\s*\('
    r'\s*Model&\s+model\s*,\s*'
    r'quadra::ParameterVector&\s+params\s*,\s*'
    r'quadra::LaplaceOptions&\s+opts\s*,\s*'
    r'const\s+quadra::OptResult&\s+fit\s*\)\s*'
    r'\{[\s\S]*?\n\}',
    re.S
)

replacement = r'''template <class Model>
Eigen::SparseMatrix<double> compute_final_random_effect_hessian(
    Model& model,
    quadra::ParameterVector& params,
    quadra::LaplaceOptions& /*opts*/,
    const quadra::OptResult& fit)
{
  // QUADRA_OPAKAPAKA_HUU_ADSCOPE_REPAIR_V1
  //
  // LaplaceResult currently stores value/gradients only. For conditional
  // random-effect SEs, rebuild the fitted AD vector, evaluate the model,
  // propagate adjoints, discover the sparse Hessian pattern, and extract H_uu
  // using Quadra's sparse Hessian extraction API.

  const std::size_t n_fixed = fit.par.size();
  const std::size_t n_random = fit.u_hat.size();
  const std::size_t n_total = n_fixed + n_random;

  std::vector<int> random_idx;
  random_idx.reserve(n_random);
  for (std::size_t i = 0; i < n_random; ++i) {
    random_idx.push_back(static_cast<int>(n_fixed + i));
  }

  quadra::ADScope scope;
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

  h_uu.makeCompressed();
  return h_uu;
}'''

s2, n = pat.subn(replacement, s, count=1)
if n != 1:
    raise SystemExit("ERROR: could not replace compute_final_random_effect_hessian function")

p.write_text(s2)
print("Patched compute_final_random_effect_hessian with ADScope + AD vector")
PY

cat > inspect_opakapaka_random_effect_selected_inverse_v3.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
exe="build/examples/opakapaka_projection"
ad_impl="build/examples/opakapaka_adgraph_global.cpp"

echo "== H_uu ADScope repair markers =="
grep -n "HUU_ADSCOPE_REPAIR\\|discover_sparse_hessian_pattern\\|extract_sparse_hessian(scope\\|scope.variable\\|model(p_full)" "$cpp" || true

if grep -q "res.H_uu" "$cpp"; then
  echo "ERROR: res.H_uu still present" >&2
  exit 1
fi

echo
echo "== Relevant core APIs =="
grep -n "discover_sparse_hessian_pattern\\|extract_sparse_hessian" core/laplace.hpp | head -40 || true

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
chmod +x inspect_opakapaka_random_effect_selected_inverse_v3.sh

echo
echo "Backups saved with suffix: huu_adscope_repair_${stamp}.bak"
echo "Run:"
echo "  ./inspect_opakapaka_random_effect_selected_inverse_v3.sh"
