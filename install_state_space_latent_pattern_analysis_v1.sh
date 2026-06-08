#!/usr/bin/env bash
set -euo pipefail

# install_state_space_latent_pattern_analysis_v1.sh
#
# Adds a focused pattern-analysis diagnostic for the state-space latent runtime.
#
# Creates:
#   examples/state_space_surplus_production/analyze_state_space_latent_pattern.cpp
#   run_state_space_surplus_latent_pattern_analysis.sh
#
# Purpose:
#   Inspect Hxx once at xhat and report:
#     - rows/cols/nnz
#     - diagonal/offdiagonal nnz
#     - max bandwidth
#     - fill ratio
#     - symmetry diagnostics
#     - detected structure
#     - backend recommendation
#     - backend name
#     - logdet agreement between sparse-H backend and direct tridiagonal values
#
# This keeps pattern analysis separate from performance timing.

src="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"
dst="examples/state_space_surplus_production/analyze_state_space_latent_pattern.cpp"

if [[ ! -f "$src" ]]; then
  echo "ERROR: missing $src"
  echo "Run the Phase 2/3 runtime installer first."
  exit 1
fi

cp "$src" "$dst"

python3 - "$dst" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

main_marker = 'int main(int argc, char** argv) {'
main_start = s.find(main_marker)
if main_start < 0:
    raise SystemExit('Could not find main function')

new_main = r'''int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();

  const Eigen::VectorXd xhat = optimize_x(data, par);
  const double joint = joint_x(data, par, xhat);
  const double grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::SparseMatrix<double> H =
      fd_tridiagonal_hessian_xx(data, par, xhat);

  quadra::laplace::StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = true;
  opts.dense_size_cutoff = 16;
  opts.banded_width_cutoff = 64;
  opts.dense_fill_ratio = 0.25;

  quadra::laplace::StructureDetector detector(opts);
  const quadra::laplace::BackendRecommendation rec = detector.Analyze(H);

  const quadra::laplace::StructureInfo info =
      quadra::laplace::InspectHessianStructure(H, opts.structure_options);

  auto backend = quadra::laplace::CreateLaplaceBackend(rec);
  backend->analyze_pattern(H);
  backend->factorize(H);

  if (!backend->is_spd()) {
    throw std::runtime_error("recommended backend reported non-SPD Hessian");
  }

  const double backend_logdet = backend->logdet();
  const double sparse_logdet = sparse_logdet_ldlt(H);

  TridiagonalValues tv = fd_tridiagonal_values_xx(data, par, xhat);
  const double value_logdet = logdet_tridiagonal_values_ldlt(tv);

  Eigen::SparseMatrix<double> H_from_values(tv.diag.size(), tv.diag.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(3 * tv.diag.size()));
  for (int i = 0; i < tv.diag.size(); ++i) {
    triplets.emplace_back(i, i, tv.diag[i]);
    if (i > 0) {
      triplets.emplace_back(i, i - 1, tv.offdiag[i - 1]);
      triplets.emplace_back(i - 1, i, tv.offdiag[i - 1]);
    }
  }
  H_from_values.setFromTriplets(triplets.begin(), triplets.end());

  const quadra::laplace::StructureInfo value_info =
      quadra::laplace::InspectHessianStructure(
          H_from_values, opts.structure_options);

  std::cout << std::fixed << std::setprecision(12);
  std::cout << "State-space latent Hessian pattern analysis\n";
  std::cout << "===========================================\n";
  std::cout << "joint = " << joint << "\n";
  std::cout << "grad_norm = " << grad_norm << "\n";
  std::cout << "x_size = " << xhat.size() << "\n";
  std::cout << "\n";

  std::cout << "Sparse H inspection\n";
  std::cout << "  rows = " << info.rows << "\n";
  std::cout << "  cols = " << info.cols << "\n";
  std::cout << "  nnz = " << info.nnz << "\n";
  std::cout << "  diagonal_nnz = " << info.diagonal_nnz << "\n";
  std::cout << "  offdiagonal_nnz = " << info.offdiagonal_nnz << "\n";
  std::cout << "  max_row_nnz = " << info.max_row_nnz << "\n";
  std::cout << "  max_bandwidth = " << info.max_bandwidth << "\n";
  std::cout << "  fill_ratio = " << info.fill_ratio << "\n";
  std::cout << "  structurally_symmetric = " << info.structurally_symmetric << "\n";
  std::cout << "  numerically_symmetric = " << info.numerically_symmetric << "\n";
  std::cout << "  max_abs_asymmetry = " << info.max_abs_asymmetry << "\n";
  std::cout << "  detected_structure = "
            << static_cast<int>(info.detected) << "\n";
  std::cout << "\n";

  std::cout << "Backend recommendation\n";
  std::cout << "  backend = " << quadra::laplace::ToString(rec.backend) << "\n";
  std::cout << "  reason = " << rec.reason << "\n";
  std::cout << "  bandwidth = " << rec.bandwidth << "\n";
  std::cout << "  fill_ratio = " << rec.fill_ratio << "\n";
  std::cout << "  supports_symbolic_reuse = "
            << rec.supports_symbolic_reuse << "\n";
  std::cout << "  backend_instance = " << backend->name() << "\n";
  std::cout << "\n";

  std::cout << "Logdet agreement\n";
  std::cout << "  sparse_logdet = " << sparse_logdet << "\n";
  std::cout << "  backend_logdet = " << backend_logdet << "\n";
  std::cout << "  value_logdet = " << value_logdet << "\n";
  std::cout << "  backend_minus_sparse = "
            << (backend_logdet - sparse_logdet) << "\n";
  std::cout << "  values_minus_sparse = "
            << (value_logdet - sparse_logdet) << "\n";
  std::cout << "\n";

  std::cout << "Value-derived H inspection\n";
  std::cout << "  nnz = " << value_info.nnz << "\n";
  std::cout << "  max_bandwidth = " << value_info.max_bandwidth << "\n";
  std::cout << "  fill_ratio = " << value_info.fill_ratio << "\n";
  std::cout << "  detected_structure = "
            << static_cast<int>(value_info.detected) << "\n";

  return 0;
}
'''

s = s[:main_start] + new_main
p.write_text(s)
PY

cat > run_state_space_surplus_latent_pattern_analysis.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -I. \
  examples/state_space_surplus_production/analyze_state_space_latent_pattern.cpp \
  -o build/examples/analyze_state_space_latent_pattern

./build/examples/analyze_state_space_latent_pattern
EOF

chmod +x run_state_space_surplus_latent_pattern_analysis.sh

cat <<'EOF'

Installed state-space latent pattern-analysis diagnostic.

Run:
  ./run_state_space_surplus_latent_pattern_analysis.sh

Expected:
  backend = tridiagonal
  max_bandwidth = 1
  logdet differences near machine precision

EOF
