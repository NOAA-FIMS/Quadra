#!/usr/bin/env bash
set -euo pipefail

# clean_state_space_latent_tridiagonal_backend_path_v1.sh
#
# Removes temporary backend diagnostics from:
#   examples/state_space_surplus_production/laplace_state_space_surplus_latent_tridiagonal.cpp
#
# Keeps the real integration:
#   H -> CreateLaplaceBackendForHessian -> backend->factorize(H) -> backend->logdet()
#
# Expected:
#   backend auto-selects tridiagonal
#   objective/logdet unchanged
#   no diagnostic print spam

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_tridiagonal.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_tridiagonal.cpp.clean_backend.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

start_marker = "out.nnz = static_cast<int>(H.nonZeros());"
end_marker = "  const double n_x = static_cast<double>(xhat.size());"

start = s.find(start_marker)
if start < 0:
    raise SystemExit("Could not find nnz/logdet block start")

end = s.find(end_marker, start)
if end < 0:
    raise SystemExit("Could not find nnz/logdet block end")

clean = '''out.nnz = static_cast<int>(H.nonZeros());

  static std::unique_ptr<quadra::laplace::LaplaceBackend> backend;
  static quadra::laplace::BackendRecommendation recommendation;

  if (!backend) {
    backend =
        quadra::laplace::CreateLaplaceBackendForHessian(
            H,
            &recommendation);
  }

  backend->factorize(H);

  if (!backend->is_spd()) {
    throw std::runtime_error(
        "Laplace backend reported non-SPD Hessian");
  }

  out.logdet = backend->logdet();

'''

s = s[:start] + clean + s[end:]
p.write_text(s)
PY

cat <<'EOF'

Cleaned state-space latent tridiagonal backend path.

Rebuild:
  c++ -std=c++17 -O2 -DNDEBUG -g     -Iexternal/Eigen     -Iexternal/LBFGSpp/include     examples/state_space_surplus_production/laplace_state_space_surplus_latent_tridiagonal.cpp     -o build/examples/laplace_state_space_surplus_latent_tridiagonal

Run:
  ./build/examples/laplace_state_space_surplus_latent_tridiagonal 20

Expected objective:
  -10.642176

EOF
