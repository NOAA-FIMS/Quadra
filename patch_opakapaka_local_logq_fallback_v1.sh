#!/usr/bin/env bash
set -euo pipefail

echo "== Patch Opakapaka: local 1D log_q fallback after LBFGS line-search stall =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
[[ -f "$cpp" ]] || { echo "ERROR: missing $cpp" >&2; exit 1; }

cp "$cpp" ".quadra_patch_backups/opakapaka_projection.cpp.local_logq_fallback_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/opakapaka_projection/opakapaka_projection.cpp")
s = p.read_text()

if "QUADRA_OPAKAPAKA_LOCAL_LOGQ_FALLBACK_V1" in s:
    print("Local log_q fallback already installed.")
    raise SystemExit(0)

# Add stdexcept/string if needed.
for inc in ["#include <stdexcept>\n", "#include <string>\n"]:
    if inc not in s and "#include <fstream>\n" in s:
        s = s.replace("#include <fstream>\n", "#include <fstream>\n" + inc, 1)

helper = r'''
// QUADRA_OPAKAPAKA_LOCAL_LOGQ_FALLBACK_V1
template <class Model>
quadra::OptResult fit_log_q_fd_newton_fallback(
    Model& model,
    quadra::ParameterVector& params,
    quadra::LaplaceOptions& opts,
    double initial_log_q)
{
  const std::vector<int> fixed_idx = {0};
  std::vector<int> random_idx;
  for (std::size_t i = 1; i < params.size(); ++i) {
    random_idx.push_back(static_cast<int>(i));
  }

  struct Eval {
    double value = std::numeric_limits<double>::infinity();
    std::vector<double> u_hat;
  };

  auto eval_at = [&](double theta) -> Eval {
    auto tmp = params;
    tmp.params.at(0).value = theta;

    Eigen::VectorXd x(1);
    x[0] = theta;

    had::ADGraph graph;
    Eval out;
    out.u_hat = quadra::solve_random_effects_laplace(
        model, tmp, x, fixed_idx, random_idx, graph);

    auto res = quadra::laplace_eval_at_u_star(
        model, tmp, fixed_idx, random_idx, x, out.u_hat, graph, opts);

    out.value = res.value;
    return out;
  };

  double theta = initial_log_q;
  Eval cur = eval_at(theta);
  double grad = std::numeric_limits<double>::infinity();
  double curv = std::numeric_limits<double>::quiet_NaN();
  int iter = 0;

  for (; iter < 25; ++iter) {
    const double h = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta)));
    const Eval left = eval_at(theta - h);
    const Eval right = eval_at(theta + h);

    if (!std::isfinite(left.value) || !std::isfinite(right.value) ||
        !std::isfinite(cur.value)) {
      break;
    }

    grad = (right.value - left.value) / (2.0 * h);
    curv = (right.value - 2.0 * cur.value + left.value) / (h * h);

    if (std::abs(grad) < 1.0e-4) {
      break;
    }
    if (!std::isfinite(curv) || curv <= 0.0) {
      break;
    }

    double step = -grad / curv;
    step = std::max(-1.0, std::min(1.0, step));

    bool accepted = false;
    for (int bt = 0; bt < 20; ++bt) {
      const double trial_theta = theta + step;
      Eval trial = eval_at(trial_theta);
      if (std::isfinite(trial.value) && trial.value <= cur.value) {
        theta = trial_theta;
        cur = std::move(trial);
        accepted = true;
        break;
      }
      step *= 0.5;
    }

    if (!accepted || std::abs(step) < 1.0e-10) {
      break;
    }
  }

  // One final centered derivative at the returned point.
  {
    const double h = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta)));
    const Eval left = eval_at(theta - h);
    const Eval right = eval_at(theta + h);
    if (std::isfinite(left.value) && std::isfinite(right.value)) {
      grad = (right.value - left.value) / (2.0 * h);
    }
  }

  params.params.at(0).value = theta;

  quadra::OptResult out;
  out.par = std::vector<double>{theta};
  out.value = cur.value;
  out.grad_norm = std::abs(grad);
  out.converged = std::abs(grad) < 1.0e-4;
  out.iterations = iter;
  out.message =
      out.converged
          ? "accepted local safeguarded one-dimensional log_q fallback after LBFGS line-search stall"
          : "local safeguarded one-dimensional log_q fallback stopped before requested tolerance";
  out.u_hat = cur.u_hat;
  return out;
}

'''

# Insert helper before int main.
m = re.search(r'(?m)^int\s+main\s*\(', s)
if not m:
    raise SystemExit("ERROR: could not find int main(")
s = s[:m.start()] + helper + "\n" + s[m.start():]

# Replace the optimize call with a guarded call. Match common formatting.
patterns = [
    r'auto\s+fit\s*=\s*quadra::optimize_lbfgs\s*\(\s*model\s*,\s*params\s*,\s*opts\s*\)\s*;',
    r'quadra::OptResult\s+fit\s*=\s*quadra::optimize_lbfgs\s*\(\s*model\s*,\s*params\s*,\s*opts\s*\)\s*;',
]

replacement = r'''quadra::OptResult fit;
  try {
    fit = quadra::optimize_lbfgs(model, params, opts);
  } catch (const std::runtime_error& e) {
    const std::string msg = e.what();
    if (msg.find("line search") == std::string::npos &&
        msg.find("sufficiently decrease") == std::string::npos) {
      throw;
    }

    std::cout << "L-BFGS line-search stall detected in Opakapaka example. "
              << "Using local safeguarded one-dimensional log_q fallback.\n";

    fit = fit_log_q_fd_newton_fallback(model, params, opts,
                                       params.params.at(0).value);
  }'''

n = 0
for pat in patterns:
    s2, n = re.subn(pat, replacement, s, count=1, flags=re.S)
    if n:
      s = s2
      break

if not n:
    raise SystemExit("ERROR: could not find quadra::optimize_lbfgs(model, params, opts) fit call")

p.write_text(s)
print("Patched Opakapaka local log_q fallback")
PY

cat > inspect_opakapaka_local_logq_fallback_v1.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
exe="build/examples/opakapaka_projection"
ad_impl="build/examples/opakapaka_adgraph_global.cpp"

echo "== Local fallback markers =="
grep -n "QUADRA_OPAKAPAKA_LOCAL_LOGQ_FALLBACK_V1\\|fit_log_q_fd_newton_fallback\\|line-search stall detected" "$cpp"

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
echo "== Level-1 output files =="
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
echo "== random_effect_uncertainty.csv preview =="
head -10 examples/opakapaka_projection/outputs/random_effect_uncertainty.csv

echo
echo "== derived_quantities.csv preview =="
head -10 examples/opakapaka_projection/outputs/derived_quantities.csv

echo
echo "== projection_uncertainty.csv preview =="
head -10 examples/opakapaka_projection/outputs/projection_uncertainty.csv
SH
chmod +x inspect_opakapaka_local_logq_fallback_v1.sh

echo
echo "Backups saved with suffix: local_logq_fallback_${stamp}.bak"
echo "Run:"
echo "  ./inspect_opakapaka_local_logq_fallback_v1.sh"
