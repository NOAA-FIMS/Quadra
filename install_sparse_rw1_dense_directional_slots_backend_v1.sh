#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.dense_directional_slots.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/dense_directional_slots_backend.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

old = """    double avg_inserts = 0.0;
    int vertices = 0;"""
new = """    double avg_inserts = 0.0;
    double dense_slots_ms = 0.0;
    double dense_speedup_vs_reuse = 0.0;
    double dense_grad_diff = 0.0;
    int vertices = 0;"""
if old not in s:
    raise SystemExit("Could not find Row counter block")
s = s.replace(old, new, 1)

if "compute_rw1_dense_slot_exact_gradient" not in s:
    marker = "struct Row {"
    idx = s.find(marker)
    if idx < 0:
        raise SystemExit("Could not find Row marker")

    helper = r'''
struct DenseSlotExactGradientResult {
    ExactGradientResult result;
    Eigen::VectorXd trace_terms;
};

Eigen::VectorXd rw1_total_u_direction(const Eigen::VectorXd& theta,
                                      const Eigen::VectorXd& uhat,
                                      int theta_index,
                                      quadra::laplace::SparseHuuFactorization& factor) {
    return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
}

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

    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);

    Eigen::VectorXd traces = Eigen::VectorXd::Zero(K);

    for (int k = 0; k < K; ++k) {
        const int j = k % 5;

        Eigen::VectorXd theta_dir = Eigen::VectorXd::Zero(5);
        theta_dir[j] = 1.0;

        const Eigen::VectorXd u_dir =
            rw1_total_u_direction(theta, uhat, j, factor);

        double trace = 0.0;

        for (int i = 0; i < m; ++i) {
            double hdot_diag = 0.0;

            if (theta_dir[1] != 0.0) {
                hdot_diag += -2.0 * inv_sigma2 * theta_dir[1];
            }

            if (theta_dir[2] != 0.0) {
                hdot_diag += lambda0 * theta_dir[2];
            }

            hdot_diag += beta * std::exp(uhat[i]) *
                         (theta_dir[4] + u_dir[i]);

            const int rw_degree = (i > 0 ? 1 : 0) + (i + 1 < m ? 1 : 0);
            if (theta_dir[3] != 0.0) {
                hdot_diag += static_cast<double>(rw_degree) *
                             lambda_rw * theta_dir[3];
            }

            trace += selected_inverse(i, i) * hdot_diag;

            if (i > 0) {
                double hdot_sub = 0.0;
                if (theta_dir[3] != 0.0) {
                    hdot_sub += -lambda_rw * theta_dir[3];
                }

                trace += 2.0 * selected_inverse(i, i - 1) * hdot_sub;
            }
        }

        traces[k] = trace;
    }

    DenseSlotExactGradientResult out;
    out.trace_terms = traces;
    out.result.objective = joint_objective + 0.5 * logdet;
    out.result.gradient = joint_grad + 0.5 * traces;
    return out;
}

'''
    s = s[:idx] + helper + s[idx:]

old = """    ExactGradientResult last_rebuild;
    ExactGradientResult last_reuse;"""
new = """    ExactGradientResult last_rebuild;
    ExactGradientResult last_reuse;
    DenseSlotExactGradientResult last_dense;"""
s = s.replace(old, new, 1)

old = """    const auto t3 = Clock::now();

    out.rebuild_ms = reuse_only"""
new = """    const auto t3 = Clock::now();

    const auto td0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        last_dense = compute_rw1_dense_slot_exact_gradient(
            m, K, theta, uhat, selected_inverse, joint_grad,
            joint_obj, logdet, factor);
    }
    const auto td1 = Clock::now();

    out.rebuild_ms = reuse_only"""
if old not in s:
    raise SystemExit("Could not find post-reuse timing block")
s = s.replace(old, new, 1)

old = """    out.avg_inserts = insert_sum / static_cast<double>(reps);
    if (reuse_only) {"""
new = """    out.avg_inserts = insert_sum / static_cast<double>(reps);
    out.dense_slots_ms = ms_between(td0, td1) / static_cast<double>(reps);
    out.dense_speedup_vs_reuse =
        out.dense_slots_ms > 0.0 ? out.reuse_ms / out.dense_slots_ms : 0.0;
    out.dense_grad_diff =
        (last_reuse.gradient - last_dense.result.gradient).cwiseAbs().maxCoeff();

    if (reuse_only) {"""
if old not in s:
    raise SystemExit("Could not find avg_inserts assignment")
s = s.replace(old, new, 1)

old = """    if (reuse_only) {
        out.grad_diff = 0.0;
        out.obj_diff = 0.0;
    } else {"""
new = """    if (reuse_only) {
        out.grad_diff = out.dense_grad_diff;
        out.obj_diff = std::abs(last_reuse.objective - last_dense.result.objective);
    } else {"""
if old in s:
    s = s.replace(old, new, 1)

old = """              << std::setw(14) << "inserts"
              << std::setw(16) << "grad diff"
              << std::setw(16) << "obj diff"
              << "\\n";"""
new = """              << std::setw(14) << "inserts"
              << std::setw(14) << "dense ms"
              << std::setw(14) << "dense spd"
              << std::setw(16) << "dense diff"
              << std::setw(16) << "grad diff"
              << std::setw(16) << "obj diff"
              << "\\n";"""
if old not in s:
    raise SystemExit("Could not find header insert point")
s = s.replace(old, new, 1)

old = """                  << std::setw(14) << r.avg_inserts
                  << std::setw(16) << r.grad_diff
                  << std::setw(16) << r.obj_diff
                  << "\\n";"""
new = """                  << std::setw(14) << r.avg_inserts
                  << std::setw(14) << r.dense_slots_ms
                  << std::setw(14) << r.dense_speedup_vs_reuse
                  << std::setw(16) << r.dense_grad_diff
                  << std::setw(16) << r.grad_diff
                  << std::setw(16) << r.obj_diff
                  << "\\n";"""
if old not in s:
    raise SystemExit("Could not find print insert point")
s = s.replace(old, new, 1)

p.write_text(s)
PYEOF

python3 /tmp/dense_directional_slots_backend.py

cat <<'EOF'

Installed experimental RW1 dense directional slots backend.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10
  ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse 10 --reuse-only

Watch:
  - dense diff should be near 0
  - dense ms should be much smaller than reuse ms

EOF
