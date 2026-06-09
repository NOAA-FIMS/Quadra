#!/usr/bin/env bash
set -euo pipefail

echo "== Repair Opakapaka random-effect selected inverse: do not use LaplaceResult::H_uu =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
lap="core/laplace.hpp"

[[ -f "$cpp" ]] || { echo "ERROR: missing $cpp" >&2; exit 1; }
[[ -f "$lap" ]] || { echo "ERROR: missing $lap" >&2; exit 1; }

cp "$cpp" ".quadra_patch_backups/opakapaka_projection.cpp.random_effect_selected_inverse_huu_repair_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

cpp = Path("examples/opakapaka_projection/opakapaka_projection.cpp")
s = cpp.read_text()

if "QUADRA_OPAKAPAKA_RANDOM_EFFECT_HUU_REPAIR_V1" in s:
    print("H_uu repair already installed.")
    raise SystemExit(0)

# Discover likely Hessian extraction function names from core/laplace.hpp.
lap = Path("core/laplace.hpp").read_text()
candidates = [
    "extract_sparse_hessian",
    "build_sparse_hessian_from_slot_workspace",
    "extract_sparse_hessian_from_graph",
]
available = [name for name in candidates if name in lap]

if not available:
    raise SystemExit(
        "ERROR: could not find a known sparse Hessian extraction function in core/laplace.hpp. "
        "Run: grep -n \"SparseMatrix<double>\\|extract_sparse_hessian\\|build_sparse_hessian\" core/laplace.hpp"
    )

# Prefer the canonical extraction path if available.
extract_name = "extract_sparse_hessian" if "extract_sparse_hessian" in available else available[0]

# Replace the bad return res.H_uu block inside compute_final_random_effect_hessian.
old = r'''  had::ADGraph graph;
  auto res = quadra::laplace_eval_at_u_star(
      model, tmp, fixed_idx, random_idx, x, fit.u_hat, graph, opts);

  return res.H_uu;
'''
new = f'''  // QUADRA_OPAKAPAKA_RANDOM_EFFECT_HUU_REPAIR_V1
  // LaplaceResult currently stores value/gradients only, so rebuild the graph
  // at the final fitted state and extract the random-effect Hessian directly.
  had::ADGraph graph;
  had::g_ADGraph = &graph;

  // Set fixed and random parameters to the fitted state.
  tmp.params.at(0).value = fit.par.at(0);
  for (std::size_t j = 0; j < fit.u_hat.size(); ++j) {{
    tmp.params.at(j + 1).value = fit.u_hat[j];
  }}

  auto objective = model(tmp);
  graph.Forward();
  had::PropagateAdjoint(objective);

  auto full_hessian = quadra::{extract_name}(
      static_cast<int>(tmp.size()));

  Eigen::SparseMatrix<double> h_uu(
      static_cast<int>(random_idx.size()),
      static_cast<int>(random_idx.size()));

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(full_hessian.nonZeros()));

  for (int k = 0; k < full_hessian.outerSize(); ++k) {{
    for (Eigen::SparseMatrix<double>::InnerIterator it(full_hessian, k); it; ++it) {{
      const int r = static_cast<int>(it.row());
      const int c = static_cast<int>(it.col());
      if (r <= 0 || c <= 0) continue;
      const int rr = r - 1;
      const int cc = c - 1;
      if (rr >= 0 && rr < static_cast<int>(random_idx.size()) &&
          cc >= 0 && cc < static_cast<int>(random_idx.size())) {{
        triplets.emplace_back(rr, cc, it.value());
      }}
    }}
  }}

  h_uu.setFromTriplets(triplets.begin(), triplets.end());
  h_uu.makeCompressed();
  return h_uu;
'''

if old not in s:
    # More flexible replacement: from "had::ADGraph graph;" through "return res.H_uu;"
    pat = re.compile(
        r'\s*had::ADGraph\s+graph\s*;\s*'
        r'auto\s+res\s*=\s*quadra::laplace_eval_at_u_star\s*\([\s\S]*?\)\s*;\s*'
        r'return\s+res\.H_uu\s*;\s*',
        re.S
    )
    s2, n = pat.subn("\n" + new, s, count=1)
    if n != 1:
        raise SystemExit("ERROR: could not find bad res.H_uu block to replace")
    s = s2
else:
    s = s.replace(old, new, 1)

cpp.write_text(s)
print(f"Patched compute_final_random_effect_hessian using quadra::{extract_name}")
PY

cat > inspect_opakapaka_random_effect_selected_inverse_v2.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
exe="build/examples/opakapaka_projection"
ad_impl="build/examples/opakapaka_adgraph_global.cpp"

echo "== H_uu repair markers =="
grep -n "RANDOM_EFFECT_HUU_REPAIR\\|extract_sparse_hessian\\|build_sparse_hessian\\|return h_uu\\|res.H_uu" "$cpp" || true

if grep -q "res.H_uu" "$cpp"; then
  echo "ERROR: res.H_uu still present" >&2
  exit 1
fi

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
chmod +x inspect_opakapaka_random_effect_selected_inverse_v2.sh

echo
echo "Backups saved with suffix: random_effect_selected_inverse_huu_repair_${stamp}.bak"
echo "Run:"
echo "  ./inspect_opakapaka_random_effect_selected_inverse_v2.sh"
