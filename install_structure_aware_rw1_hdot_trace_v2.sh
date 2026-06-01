#!/usr/bin/env bash
set -euo pipefail

# install_structure_aware_rw1_hdot_trace_v2.sh
#
# Promotes the validated dense RW1 Hdot/trace backend into:
#   core/laplace/structure_aware_rw1_hdot.hpp
#
# Then patches the benchmark to call that reusable helper.

mkdir -p core/laplace .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.structure_aware_rw1_v2.$(date +%Y%m%d_%H%M%S).bak"

cat > core/laplace/structure_aware_rw1_hdot.hpp <<'EOF'
#ifndef QUADRA_LAPLACE_STRUCTURE_AWARE_RW1_HDOT_HPP
#define QUADRA_LAPLACE_STRUCTURE_AWARE_RW1_HDOT_HPP

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>
#include <utility>

#include "sparse_huu_factorization.hpp"

namespace quadra {
namespace laplace {

struct StructureAwareRw1HdotTraceResult {
  Eigen::VectorXd trace_terms;
  Eigen::VectorXd gradient;
  double objective = 0.0;
};

template <class FuThetaColumnProvider>
Eigen::VectorXd rw1_total_u_direction(
    const Eigen::VectorXd &theta,
    const Eigen::VectorXd &uhat,
    const int theta_index,
    SparseHuuFactorization &factor,
    FuThetaColumnProvider &&f_u_theta_column) {
  (void)theta;

  if (theta_index < 0) {
    throw std::invalid_argument("rw1_total_u_direction: theta_index < 0");
  }

  return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
}

template <class SelectedInverseAccessor, class FuThetaColumnProvider>
Eigen::VectorXd rw1_structure_aware_hdot_traces(
    const int m,
    const int n_directions,
    const Eigen::VectorXd &theta,
    const Eigen::VectorXd &uhat,
    SelectedInverseAccessor &&selected_inverse,
    SparseHuuFactorization &factor,
    FuThetaColumnProvider &&f_u_theta_column) {
  if (m < 0) {
    throw std::invalid_argument("rw1_structure_aware_hdot_traces: m < 0");
  }
  if (n_directions < 0) {
    throw std::invalid_argument(
        "rw1_structure_aware_hdot_traces: n_directions < 0");
  }
  if (theta.size() < 5) {
    throw std::invalid_argument(
        "rw1_structure_aware_hdot_traces: theta must have at least 5 values");
  }
  if (uhat.size() != m) {
    throw std::invalid_argument(
        "rw1_structure_aware_hdot_traces: uhat size does not match m");
  }

  const double lambda0 = std::exp(theta[2]);
  const double beta = std::exp(theta[4]);

  Eigen::VectorXd traces = Eigen::VectorXd::Zero(n_directions);

  for (int k = 0; k < n_directions; ++k) {
    const int theta_index = k % 5;

    Eigen::VectorXd theta_direction = Eigen::VectorXd::Zero(5);
    theta_direction[theta_index] = 1.0;

    const Eigen::VectorXd u_direction =
        rw1_total_u_direction(theta, uhat, theta_index, factor,
                              std::forward<FuThetaColumnProvider>(
                                  f_u_theta_column));

    double trace = 0.0;

    for (int i = 0; i < m; ++i) {
      double hdot_diag = 0.0;

      // Direct log_lambda0 contribution.
      if (theta_direction[2] != 0.0) {
        hdot_diag += lambda0 * theta_direction[2];
      }

      // Total direction contribution through beta * exp(u_i).
      hdot_diag += beta * std::exp(uhat[i]) *
                   (theta_direction[4] + u_direction[i]);

      trace += selected_inverse(i, i) * hdot_diag;

      if (i > 0) {
        // Current validated HAD Hdot convention has no direct RW off-diagonal
        // term in this dense-slot comparison.
        const double hdot_subdiag = 0.0;
        trace += 2.0 * selected_inverse(i, i - 1) * hdot_subdiag;
      }
    }

    traces[k] = trace;
  }

  return traces;
}

template <class SelectedInverseAccessor, class FuThetaColumnProvider>
StructureAwareRw1HdotTraceResult rw1_structure_aware_exact_gradient(
    const int m,
    const int n_directions,
    const Eigen::VectorXd &theta,
    const Eigen::VectorXd &uhat,
    SelectedInverseAccessor &&selected_inverse,
    const Eigen::VectorXd &joint_gradient,
    const double joint_objective,
    const double logdet_huu,
    SparseHuuFactorization &factor,
    FuThetaColumnProvider &&f_u_theta_column) {
  StructureAwareRw1HdotTraceResult out;
  out.trace_terms = rw1_structure_aware_hdot_traces(
      m, n_directions, theta, uhat,
      std::forward<SelectedInverseAccessor>(selected_inverse),
      factor,
      std::forward<FuThetaColumnProvider>(f_u_theta_column));

  out.objective = joint_objective + 0.5 * logdet_huu;
  out.gradient = joint_gradient + 0.5 * out.trace_terms;
  return out;
}

}  // namespace laplace
}  // namespace quadra

#endif  // QUADRA_LAPLACE_STRUCTURE_AWARE_RW1_HDOT_HPP
EOF

cat > /tmp/patch_benchmark_structure_aware_rw1_v2.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

include_line = '#include "../core/laplace/structure_aware_rw1_hdot.hpp"'
if include_line not in s:
    s = s.replace(
        '#include "../core/laplace/sparse_huu_factorization.hpp"',
        '#include "../core/laplace/sparse_huu_factorization.hpp"\n' + include_line,
        1,
    )

start = s.find("DenseSlotExactGradientResult compute_rw1_dense_slot_exact_gradient(")
if start < 0:
    raise SystemExit("dense helper not found")

brace = s.find("{", start)
depth = 0
end = None
for i in range(brace, len(s)):
    if s[i] == "{":
        depth += 1
    elif s[i] == "}":
        depth -= 1
        if depth == 0:
            end = i + 1
            break

if end is None:
    raise SystemExit("could not find end of dense helper")

new_helper = """
DenseSlotExactGradientResult compute_rw1_dense_slot_exact_gradient(
    int m,
    int K,
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& uhat,
    const SelectedTridiagonalInverse& selected_inverse,
    const Eigen::VectorXd& joint_grad,
    double joint_objective,
    double logdet,
    quadra::laplace::SparseHuuFactorization& factor) {

    const auto structure_aware =
        quadra::laplace::rw1_structure_aware_exact_gradient(
            m,
            K,
            theta,
            uhat,
            selected_inverse,
            joint_grad,
            joint_objective,
            logdet,
            factor,
            [](const Eigen::VectorXd& theta_arg,
               const Eigen::VectorXd& uhat_arg,
               int theta_index) {
                return f_u_theta_column(theta_arg, uhat_arg, theta_index);
            });

    DenseSlotExactGradientResult out;
    out.trace_terms = structure_aware.trace_terms;
    out.result.objective = structure_aware.objective;
    out.result.gradient = structure_aware.gradient;
    return out;
}
"""

s = s[:start] + new_helper + s[end:]
p.write_text(s)
PYEOF

python3 /tmp/patch_benchmark_structure_aware_rw1_v2.py

cat <<'EOF'

Installed structure-aware RW1 Hdot/trace helper v2.

Added:
  core/laplace/structure_aware_rw1_hdot.hpp

Patched:
  benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 1
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Expected:
  dense diff = 0
  grad diff = 0
  dense ms remains tiny

EOF
