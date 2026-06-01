#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_reusable_had_graph_workspace.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

if [[ ! -f core/had_graph_workspace.hpp ]]; then
  echo "ERROR: missing core/had_graph_workspace.hpp"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_reusable_had_graph_workspace.cpp.use_had_graph_workspace.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_patch_reusable_rw1_workspace.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_reusable_had_graph_workspace.cpp")
s = p.read_text()

if '#include "../core/had_graph_workspace.hpp"' not in s:
    s = s.replace(
        '#include "../core/had_quadra.hpp"\n',
        '#include "../core/had_quadra.hpp"\n#include "../core/had_graph_workspace.hpp"\n'
    )

start = s.find("struct ReusableRw1HadWorkspace {")
if start < 0:
    raise SystemExit("Could not find local ReusableRw1HadWorkspace struct")

brace = s.find("{", start)
depth = 0
end = None
for i in range(brace, len(s)):
    if s[i] == "{":
        depth += 1
    elif s[i] == "}":
        depth -= 1
        if depth == 0:
            j = i + 1
            while j < len(s) and s[j].isspace():
                j += 1
            end = j + 1 if j < len(s) and s[j] == ";" else i + 1
            break

if end is None:
    raise SystemExit("Could not find end of ReusableRw1HadWorkspace struct")

replacement = """
struct ReusableRw1HadWorkspace {
    int m = 0;
    int K = 0;
    quadra::HadGraphWorkspace had_workspace;
    std::vector<AReal> x;
    AReal f;

    ReusableRw1HadWorkspace(int m_, int K_,
                            const Eigen::VectorXd& theta,
                            const Eigen::VectorXd& uhat)
        : m(m_), K(K_), x(static_cast<size_t>(5 + m_)) {
        build(theta, uhat);
    }

    void build(const Eigen::VectorXd& theta,
               const Eigen::VectorXd& uhat) {
        f = had_workspace.Build([&]() {
            for (int j = 0; j < 5; ++j) {
                x[static_cast<size_t>(j)] = AReal(theta[j]);
            }
            for (int i = 0; i < m; ++i) {
                x[static_cast<size_t>(5 + i)] = AReal(uhat[i]);
            }

            SparseRw1Objective objective{m};
            return objective(x);
        });

        had_workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));
    }

    void propagate_base_adjoint() {
        had_workspace.PropagateAdjoint(f.varId);
    }

    void seed_directions(const Eigen::VectorXd& theta,
                         const Eigen::VectorXd& uhat) {
        const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
        quadra::laplace::SparseHuuFactorization factor(H);

        auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
            return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
        };

        had_workspace.Activate();
        had_workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));

        for (int k = 0; k < K; ++k) {
            const int theta_index = k % 5;
            const Eigen::VectorXd udir = direction_provider(theta_index);

            for (int j = 0; j < 5; ++j) {
                had::SetARealDotBatch(x[static_cast<size_t>(j)], k,
                                      j == theta_index ? 1.0 : 0.0);
            }
            for (int i = 0; i < m; ++i) {
                had::SetARealDotBatch(x[static_cast<size_t>(5 + i)], k,
                                      udir[i]);
            }
        }
    }

    void propagate_directional_batch() {
        had_workspace.PropagateAdjointDirectionalBatch();
    }

    double extract_checksum() {
        had_workspace.Activate();

        double checksum = 0.0;

        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < m; ++i) {
                checksum += had::GetAdjointDotBatch(
                    x[static_cast<size_t>(5 + i)],
                    x[static_cast<size_t>(5 + i)],
                    k);

                if (i > 0) {
                    checksum += had::GetAdjointDotBatch(
                        x[static_cast<size_t>(5 + i)],
                        x[static_cast<size_t>(5 + i - 1)],
                        k);
                }
            }
        }

        return checksum;
    }

    int vertex_count() const {
        return static_cast<int>(had_workspace.Graph().vertices.size());
    }
};
"""

s = s[:start] + replacement + s[end:]

s = s.replace(
    "had::PropagateAdjointDirectionalBatch();\n        reuse_checksum += workspace.extract_checksum();",
    "workspace.propagate_directional_batch();\n        reuse_checksum += workspace.extract_checksum();"
)

s = s.replace(
    "out.vertices = static_cast<int>(workspace.graph.vertices.size());",
    "out.vertices = workspace.vertex_count();"
)

p.write_text(s)
PYEOF

python3 /tmp/quadra_patch_reusable_rw1_workspace.py

cat <<'EOF'

Refactored reusable RW1 benchmark to use quadra::HadGraphWorkspace.

Run:
  ./run_had_graph_workspace_test.sh
  ./run_sparse_rw1_reusable_had_graph_workspace.sh 10

EOF
