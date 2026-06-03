#!/usr/bin/env bash
set -euo pipefail

target="examples/age_structured_recruitment/benchmark_age_structured_no_plus_analytic_banded.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_age_structured_no_plus_analytic_banded_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/benchmark_age_structured_no_plus_analytic_banded.cpp.newton.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("examples/age_structured_recruitment/benchmark_age_structured_no_plus_analytic_banded.cpp")
s = p.read_text()

start = s.find("Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {")
if start < 0:
    raise SystemExit("Could not find optimize_x start")

marker = "\ndouble sparse_logdet_ldlt"
end = s.find(marker, start)
if end < 0:
    raise SystemExit("Could not find optimize_x end marker")

new = r'''Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(data.n_years);

  double f = eval_all(data, par, x).objective;

  for (int iter = 0; iter < 50; ++iter) {
    const EvalAll e = eval_all(data, par, x);
    f = e.objective;

    const double gnorm = e.gradient.norm();
    if (gnorm < 1e-8) {
      break;
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(e.hessian);

    if (ldlt.info() != Eigen::Success) {
      throw std::runtime_error("Newton Hessian factorization failed");
    }

    const Eigen::VectorXd step_direction = ldlt.solve(e.gradient);

    if (ldlt.info() != Eigen::Success || !step_direction.allFinite()) {
      throw std::runtime_error("Newton Hessian solve failed");
    }

    double step = 1.0;
    bool accepted = false;

    for (int ls = 0; ls < 30; ++ls) {
      const Eigen::VectorXd candidate = x - step * step_direction;
      const double f_candidate = eval_all(data, par, candidate).objective;

      if (std::isfinite(f_candidate) && f_candidate < f) {
        x = candidate;
        f = f_candidate;
        accepted = true;
        break;
      }

      step *= 0.5;
    }

    if (!accepted) {
      if (gnorm < 1e-5) {
        break;
      }

      throw std::runtime_error("Newton line search failed");
    }
  }

  return x;
}

'''

s = s[:start] + new + s[end+1:]
p.write_text(s)
PYEOF

cat <<'EOF'

Replaced LBFGS++ random-effects optimizer with damped sparse Newton.

Run:
  ./run_quadra_analytic_vs_tmb_age_structured_no_plus_benchmark.sh 10 25,50,100,250,500,1000 10

Quick:
  ./run_quadra_analytic_vs_tmb_age_structured_no_plus_benchmark.sh 3 25,50 10

EOF
