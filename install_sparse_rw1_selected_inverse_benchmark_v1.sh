#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_selected_inverse_benchmark_v1.sh
#
# Refactors the sparse RW1 exact-gradient graph reuse benchmark so the reuse
# path uses selected inverse entries instead of dense Hinv access.
#
# For this benchmark, we compute only the tridiagonal selected inverse entries:
#   H^{-1}_{i,i}
#   H^{-1}_{i,i-1}
#
# This is still implemented by sparse solves against unit vectors for benchmark
# simplicity, but it avoids materializing/using a dense Hinv in the reuse path
# and validates the TraceTermsSelectedInverse API in a nonzero path.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.selected_inverse.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_selected_inverse_benchmark_patch.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

if "struct SelectedTridiagonalInverse" not in s:
    marker = "struct Row {"
    idx = s.find(marker)
    if idx < 0:
        raise SystemExit("Could not find Row marker")

    block = """
struct SelectedTridiagonalInverse {
    Eigen::VectorXd diag;
    Eigen::VectorXd subdiag;

    double operator()(int row, int col) const {
        if (row == col) {
            return diag[row];
        }

        const int hi = std::max(row, col);
        const int lo = std::min(row, col);

        if (hi == lo + 1) {
            return subdiag[hi - 1];
        }

        return 0.0;
    }
};

SelectedTridiagonalInverse compute_selected_tridiagonal_inverse(
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& ldlt,
    int m) {
    SelectedTridiagonalInverse selected;
    selected.diag = Eigen::VectorXd::Zero(m);
    selected.subdiag = Eigen::VectorXd::Zero(std::max(0, m - 1));

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m);

    for (int j = 0; j < m; ++j) {
        rhs.setZero();
        rhs[j] = 1.0;

        const Eigen::VectorXd col = ldlt.solve(rhs);

        selected.diag[j] = col[j];

        if (j > 0) {
            selected.subdiag[j - 1] = col[j - 1];
        }
    }

    return selected;
}

"""
    s = s[:idx] + block + s[idx:]

# Change adapter exact_gradient signature and body to accept selected inverse accessor.
old_sig = """    ExactGradientResult exact_gradient(const Eigen::MatrixXd& Hinv,
                                       const Eigen::VectorXd& joint_grad,
                                       double joint_objective,
                                       double logdet) {"""
new_sig = """    template <class SelectedInverseAccessor>
    ExactGradientResult exact_gradient(SelectedInverseAccessor&& selected_inverse,
                                       const Eigen::VectorXd& joint_grad,
                                       double joint_objective,
                                       double logdet) {"""
if old_sig in s:
    s = s.replace(old_sig, new_sig, 1)

old_lambda = """                [&](int row, int col) {
                    return Hinv(row, col);
                },"""
new_lambda = """                std::forward<SelectedInverseAccessor>(selected_inverse),"""
if old_lambda in s:
    s = s.replace(old_lambda, new_lambda, 1)

# Replace dense Hinv computation block with selected inverse computation.
old_hinv = """    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(H);
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
    const Eigen::MatrixXd Hinv = ldlt.solve(I);"""
new_hinv = """    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(H);
    SelectedTridiagonalInverse selected_inverse =
        compute_selected_tridiagonal_inverse(ldlt, m);

    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
    const Eigen::MatrixXd Hinv = ldlt.solve(I);  // Rebuild reference only."""
if old_hinv in s:
    s = s.replace(old_hinv, new_hinv, 1)
else:
    print("WARNING: dense Hinv block not found or already modified")

# Replace reuse exact_gradient call.
old_call = """        last_reuse = workspace.exact_gradient(
            Hinv, joint_grad, joint_obj, logdet);"""
new_call = """        last_reuse = workspace.exact_gradient(
            selected_inverse, joint_grad, joint_obj, logdet);"""
if old_call in s:
    s = s.replace(old_call, new_call, 1)
else:
    print("WARNING: reuse exact_gradient call not found")

p.write_text(s)
PYEOF

python3 /tmp/quadra_selected_inverse_benchmark_patch.py

cat <<'EOF'

Installed selected-inverse benchmark refactor.

Patched:
  benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Notes:
  - The rebuild reference still uses dense Hinv for comparison.
  - The reuse path uses the selected inverse accessor.

EOF
