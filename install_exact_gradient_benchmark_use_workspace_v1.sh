#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Install/run install_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh first."
  exit 1
fi

if [[ ! -f core/laplace/exact_gradient_workspace.hpp ]]; then
  echo "ERROR: missing core/laplace/exact_gradient_workspace.hpp"
  echo "Install/run install_exact_gradient_workspace_v1.sh first."
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.use_exact_gradient_workspace.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_patch_exact_gradient_benchmark_use_workspace.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

if '#include "../core/laplace/exact_gradient_workspace.hpp"' not in s:
    s = s.replace(
        '#include "../core/had_quadra.hpp"\n',
        '#include "../core/had_quadra.hpp"\n#include "../core/laplace/exact_gradient_workspace.hpp"\n'
    )

if "struct ExactGradientWorkspaceRw1Adapter" not in s:
    marker = "struct Row {"
    idx = s.find(marker)
    if idx < 0:
        raise SystemExit("Could not find Row marker")

    adapter = """
struct ExactGradientWorkspaceRw1Adapter {
    int m = 0;
    int K = 0;
    quadra::laplace::ExactGradientWorkspace workspace;
    std::vector<AReal> theta_vars;
    std::vector<AReal> random_vars;
    AReal f;

    ExactGradientWorkspaceRw1Adapter(int m_, int K_,
                                     const Eigen::VectorXd& theta,
                                     const Eigen::VectorXd& uhat)
        : m(m_),
          K(K_),
          theta_vars(5),
          random_vars(static_cast<std::size_t>(m_)) {
        build(theta, uhat);
    }

    void build(const Eigen::VectorXd& theta,
               const Eigen::VectorXd& uhat) {
        f = workspace.Build(
            [&]() {
                for (int j = 0; j < 5; ++j) {
                    theta_vars[static_cast<std::size_t>(j)] = AReal(theta[j]);
                }
                for (int i = 0; i < m; ++i) {
                    random_vars[static_cast<std::size_t>(i)] = AReal(uhat[i]);
                }

                std::vector<AReal> x;
                x.reserve(static_cast<std::size_t>(5 + m));

                for (auto& v : theta_vars) {
                    x.push_back(v);
                }
                for (auto& v : random_vars) {
                    x.push_back(v);
                }

                SparseRw1Objective objective{m};
                return objective(x);
            },
            &theta_vars,
            &random_vars);

        workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));
    }

    void propagate_base_adjoint() {
        workspace.PropagateBaseAdjoint();
    }

    void seed_directions(const Eigen::VectorXd& theta,
                         const Eigen::VectorXd& uhat,
                         quadra::laplace::SparseHuuFactorization& factor) {
        workspace.SeedTotalDirections(
            static_cast<std::size_t>(K),
            [&](std::size_t k,
                Eigen::VectorXd& theta_direction,
                Eigen::VectorXd& random_direction) {
                const int theta_index = static_cast<int>(k % 5);

                theta_direction = Eigen::VectorXd::Zero(5);
                theta_direction[theta_index] = 1.0;

                random_direction =
                    -factor.solve(
                        f_u_theta_column(theta, uhat, theta_index));
            });
    }

    void propagate_directional_batch() {
        workspace.PropagateDirectionalBatch();
    }

    ExactGradientResult exact_gradient(const Eigen::MatrixXd& Hinv,
                                       const Eigen::VectorXd& joint_grad,
                                       double joint_objective,
                                       double logdet) {
        const auto pattern =
            quadra::laplace::MakeTridiagonalHdotPattern(m);

        const Eigen::VectorXd traces =
            workspace.TraceTerms(Hinv, pattern);

        ExactGradientResult out;
        out.objective = joint_objective + 0.5 * logdet;
        out.gradient = joint_grad + 0.5 * traces;
        return out;
    }

    int vertex_count() const {
        return static_cast<int>(
            workspace.HadWorkspace().Graph().vertices.size());
    }
};

"""
    s = s[:idx] + adapter + s[idx:]

old = """    ReusableRw1HadWorkspace workspace(m, K, theta, uhat);
    out.vertices = workspace.vertex_count();

    const auto t2 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        workspace.propagate_base_adjoint();
        workspace.seed_directions(theta, uhat, factor);
        workspace.propagate_directional_batch();

        last_reuse = compute_exact_gradient_from_current_graph(
            workspace.x, m, K, Hinv, joint_grad, joint_obj, logdet);
    }
    const auto t3 = Clock::now();"""

new = """    ExactGradientWorkspaceRw1Adapter workspace(m, K, theta, uhat);
    out.vertices = workspace.vertex_count();

    const auto t2 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        workspace.propagate_base_adjoint();
        workspace.seed_directions(theta, uhat, factor);
        workspace.propagate_directional_batch();

        last_reuse = workspace.exact_gradient(
            Hinv, joint_grad, joint_obj, logdet);
    }
    const auto t3 = Clock::now();"""

if old not in s:
    raise SystemExit("Could not find reuse block to replace. The benchmark may already be modified.")

s = s.replace(old, new, 1)

p.write_text(s)
PYEOF

python3 /tmp/quadra_patch_exact_gradient_benchmark_use_workspace.py

cat <<'EOF'

Refactored exact-gradient graph-reuse benchmark to use ExactGradientWorkspace.

Run:
  ./run_exact_gradient_workspace_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
