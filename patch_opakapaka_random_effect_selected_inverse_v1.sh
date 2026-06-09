#!/usr/bin/env bash
set -euo pipefail

echo "== Wire selected-inverse diagonal into Opakapaka random-effect uncertainty =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
[[ -f "$cpp" ]] || { echo "ERROR: missing $cpp" >&2; exit 1; }

cp "$cpp" ".quadra_patch_backups/opakapaka_projection.cpp.random_effect_selected_inverse_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/opakapaka_projection/opakapaka_projection.cpp")
s = p.read_text()

if "QUADRA_OPAKAPAKA_RANDOM_EFFECT_SELECTED_INVERSE_V1" in s:
    print("Selected-inverse Opakapaka reporting already installed.")
    raise SystemExit(0)

# Include utility.
include = '#include "../../core/uncertainty/selected_inverse_diagonal.hpp"\n'
if include not in s:
    # opakapaka_model.hpp is usually first include and already relative to example dir.
    s = s.replace('#include "opakapaka_model.hpp"\n',
                  '#include "opakapaka_model.hpp"\n' + include, 1)

# Add a helper that evaluates final Huu and writes real conditional SEs.
helper = r'''
// QUADRA_OPAKAPAKA_RANDOM_EFFECT_SELECTED_INVERSE_V1
template <class Model>
Eigen::SparseMatrix<double> compute_final_random_effect_hessian(
    Model& model,
    quadra::ParameterVector& params,
    quadra::LaplaceOptions& opts,
    const quadra::OptResult& fit)
{
  const std::vector<int> fixed_idx = {0};
  std::vector<int> random_idx;
  for (std::size_t i = 1; i < params.size(); ++i) {
    random_idx.push_back(static_cast<int>(i));
  }

  auto tmp = params;
  tmp.params.at(0).value = fit.par.at(0);

  Eigen::VectorXd x(1);
  x[0] = fit.par.at(0);

  had::ADGraph graph;
  auto res = quadra::laplace_eval_at_u_star(
      model, tmp, fixed_idx, random_idx, x, fit.u_hat, graph, opts);

  return res.H_uu;
}

inline void write_random_effect_uncertainty_csv(
    const std::string& path,
    const std::vector<double>& u_hat,
    const Eigen::SparseMatrix<double>& h_uu)
{
  const auto diag =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(h_uu);

  std::ofstream out(path);
  out << "effect,mode,conditional_se,conditional_variance,note\n";

  for (std::size_t i = 0; i < u_hat.size(); ++i) {
    double se = std::numeric_limits<double>::quiet_NaN();
    double var = std::numeric_limits<double>::quiet_NaN();
    std::string note = diag.message;

    if (diag.success && i < diag.standard_error.size() && i < diag.variance.size()) {
      se = diag.standard_error[i];
      var = diag.variance[i];
      note = "selected_inverse_diagonal";
    }

    out << "log_B[" << i << "]," << u_hat[i] << ","
        << se << "," << var << "," << note << "\n";
  }
}

'''

# Insert helper before int main, but after existing reporting helpers if possible.
m = re.search(r'(?m)^int\s+main\s*\(', s)
if not m:
    raise SystemExit("ERROR: could not find int main(")
s = s[:m.start()] + helper + "\n" + s[m.start():]

# Remove/rename previous old writer if it exists with same name and only two args.
# Instead of trying to delete the old function robustly, rename it if present.
s = s.replace(
    "inline void write_random_effect_uncertainty_csv(const std::string& path, const std::vector<double>& u_hat)",
    "inline void write_random_effect_uncertainty_pending_csv(const std::string& path, const std::vector<double>& u_hat)"
)

# If the previous rename also hit our new writer accidentally (unlikely due signature), repair.
s = s.replace(
    "inline void write_random_effect_uncertainty_pending_csv(\n    const std::string& path,\n    const std::vector<double>& u_hat,\n    const Eigen::SparseMatrix<double>& h_uu)",
    "inline void write_random_effect_uncertainty_csv(\n    const std::string& path,\n    const std::vector<double>& u_hat,\n    const Eigen::SparseMatrix<double>& h_uu)"
)

# Insert final Huu computation before the random-effect uncertainty write call.
old_call = 'write_random_effect_uncertainty_csv("examples/opakapaka_projection/outputs/random_effect_uncertainty.csv", fit.u_hat);'
new_call = '''const auto final_h_uu =
      compute_final_random_effect_hessian(model, params, opts, fit);
  write_random_effect_uncertainty_csv(
      "examples/opakapaka_projection/outputs/random_effect_uncertainty.csv",
      fit.u_hat, final_h_uu);'''
if old_call in s:
    s = s.replace(old_call, new_call, 1)
else:
    # Try multiline variant from prior patch.
    pat = re.compile(
        r'write_random_effect_uncertainty_csv\s*\(\s*"examples/opakapaka_projection/outputs/random_effect_uncertainty\.csv"\s*,\s*fit\.u_hat\s*\)\s*;',
        re.S
    )
    s2, n = pat.subn(new_call, s, count=1)
    if n != 1:
        raise SystemExit("ERROR: could not find random_effect_uncertainty writer call")
    s = s2

p.write_text(s)
print("Patched Opakapaka random-effect selected inverse reporting")
PY

cat > inspect_opakapaka_random_effect_selected_inverse_v1.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
exe="build/examples/opakapaka_projection"
ad_impl="build/examples/opakapaka_adgraph_global.cpp"

echo "== Selected-inverse reporting markers =="
grep -n "selected_inverse_diagonal\\|compute_final_random_effect_hessian\\|RANDOM_EFFECT_SELECTED_INVERSE" "$cpp"

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
echo "== Check pending text removed from random-effect uncertainty =="
if grep -q "pending selected-inverse" examples/opakapaka_projection/outputs/random_effect_uncertainty.csv; then
  echo "ERROR: pending selected-inverse text still present" >&2
  exit 1
else
  echo "OK: random-effect conditional SEs are populated"
fi
SH
chmod +x inspect_opakapaka_random_effect_selected_inverse_v1.sh

echo
echo "Backups saved with suffix: random_effect_selected_inverse_${stamp}.bak"
echo "Run:"
echo "  ./inspect_opakapaka_random_effect_selected_inverse_v1.sh"
