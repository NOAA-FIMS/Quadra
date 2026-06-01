#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

header="core/laplace/exact_gradient_workspace.hpp"
test_file="tests/test_exact_gradient_workspace.cpp"
bench_file="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  exit 1
fi

if [[ ! -f "$test_file" ]]; then
  echo "ERROR: missing $test_file"
  exit 1
fi

cp "$header" ".quadra_patch_backups/exact_gradient_workspace.hpp.assemble_result.$(date +%Y%m%d_%H%M%S).bak"
cp "$test_file" ".quadra_patch_backups/test_exact_gradient_workspace.cpp.assemble_result.$(date +%Y%m%d_%H%M%S).bak"
if [[ -f "$bench_file" ]]; then
  cp "$bench_file" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.assemble_result.$(date +%Y%m%d_%H%M%S).bak"
fi

cat > /tmp/quadra_exact_gradient_assemble_patch.py <<'PYEOF'
from pathlib import Path

header = Path("core/laplace/exact_gradient_workspace.hpp")
s = header.read_text()

if "struct ExactGradientEvaluation" not in s:
    marker = "struct SparseHdotPatternEntry {"
    idx = s.find(marker)
    if idx < 0:
        raise SystemExit("Could not find SparseHdotPatternEntry marker")

    result_struct = """struct ExactGradientEvaluation {
    double objective = 0.0;
    Eigen::VectorXd gradient;
    Eigen::VectorXd trace_terms;
};

"""
    s = s[:idx] + result_struct + s[idx:]

if "AssembleExactGradient" not in s:
    marker = "    HadGraphWorkspace& HadWorkspace() { return had_workspace_; }"
    idx = s.find(marker)
    if idx < 0:
        raise SystemExit("Could not find HadWorkspace marker")

    method = """
    template <class SelectedInverseAccessor>
    ExactGradientEvaluation AssembleExactGradient(
        double joint_objective,
        double logdet_huu,
        const Eigen::VectorXd& joint_envelope_gradient,
        SelectedInverseAccessor&& selected_inverse,
        const std::vector<SparseHdotPatternEntry>& pattern) {
        RequireBuilt();

        if (joint_envelope_gradient.size() !=
            static_cast<int>(n_directions_)) {
            throw std::invalid_argument(
                "ExactGradientWorkspace::AssembleExactGradient gradient dimension mismatch.");
        }

        ExactGradientEvaluation out;
        out.objective = joint_objective + 0.5 * logdet_huu;
        out.trace_terms =
            TraceTermsSelectedInverse(
                std::forward<SelectedInverseAccessor>(selected_inverse),
                pattern);
        out.gradient = joint_envelope_gradient + 0.5 * out.trace_terms;

        return out;
    }

"""
    s = s[:idx] + method + s[idx:]

header.write_text(s)

test = Path("tests/test_exact_gradient_workspace.cpp")
t = test.read_text()

if "AssembleExactGradient" not in t:
    old = """    const Eigen::MatrixXd hdot0 = workspace.ExtractHdotDense(0, pattern);"""

    new = """    const Eigen::VectorXd joint_gradient = Eigen::VectorXd::Zero(2);

    const auto assembled =
        workspace.AssembleExactGradient(
            12.5,
            1.25,
            joint_gradient,
            [&](int row, int col) {
                return Hinv(row, col);
            },
            pattern);

    if (std::abs(assembled.objective - 13.125) > 1.0e-12) {
        throw std::runtime_error(
            "AssembleExactGradient returned wrong objective");
    }

    if (assembled.gradient.size() != 2 ||
        assembled.trace_terms.size() != 2) {
        throw std::runtime_error(
            "AssembleExactGradient returned wrong vector sizes");
    }

    if ((assembled.gradient - 0.5 * traces).cwiseAbs().maxCoeff() >
        1.0e-12) {
        throw std::runtime_error(
            "AssembleExactGradient gradient does not match trace assembly");
    }

    const Eigen::MatrixXd hdot0 = workspace.ExtractHdotDense(0, pattern);"""

    if old not in t:
        raise SystemExit("Could not find insertion point in test_exact_gradient_workspace.cpp")

    t = t.replace(old, new, 1)
    test.write_text(t)

bench = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
if bench.exists():
    b = bench.read_text()

    old = """    ExactGradientResult exact_gradient(const Eigen::MatrixXd& Hinv,
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
    }"""

    new = """    ExactGradientResult exact_gradient(const Eigen::MatrixXd& Hinv,
                                       const Eigen::VectorXd& joint_grad,
                                       double joint_objective,
                                       double logdet) {
        const auto pattern =
            quadra::laplace::MakeTridiagonalHdotPattern(m);

        const auto assembled =
            workspace.AssembleExactGradient(
                joint_objective,
                logdet,
                joint_grad,
                [&](int row, int col) {
                    return Hinv(row, col);
                },
                pattern);

        ExactGradientResult out;
        out.objective = assembled.objective;
        out.gradient = assembled.gradient;
        return out;
    }"""

    if old in b:
        b = b.replace(old, new, 1)
        bench.write_text(b)
PYEOF

python3 /tmp/quadra_exact_gradient_assemble_patch.py

cat <<'EOF'

Installed ExactGradientWorkspace exact-gradient assembly API.

Patched:
  core/laplace/exact_gradient_workspace.hpp
  tests/test_exact_gradient_workspace.cpp
  benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp (if matching block found)

Run:
  ./run_exact_gradient_workspace_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
