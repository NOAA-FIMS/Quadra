#!/usr/bin/env bash
set -euo pipefail

FILE="examples/state_space_surplus_production/laplace_state_space_surplus_latent_tridiagonal.cpp"

cp "$FILE" "${FILE}.backend_diag.bak"

python3 <<'PY'
from pathlib import Path

p = Path("examples/state_space_surplus_production/laplace_state_space_surplus_latent_tridiagonal.cpp")
txt = p.read_text()

if "[backend create]" in txt:
    print("Diagnostics already installed.")
    raise SystemExit(0)

needle = """out.nnz = static_cast<int>(H.nonZeros());
  out.logdet = sparse_logdet_ldlt(H);"""

replacement = """out.nnz = static_cast<int>(H.nonZeros());

  using namespace quadra::laplace;

  static std::unique_ptr<LaplaceBackend> backend;
  static BackendRecommendation recommendation;
  static int backend_creations = 0;

  if (!backend) {
    ++backend_creations;

    backend =
        CreateLaplaceBackendForHessian(
            H,
            &recommendation);

    std::cout << "[backend create] "
              << backend_creations
              << "\\n";

    std::cout << "[backend] "
              << backend->name()
              << "\\n";
  }

  const double old_logdet =
      sparse_logdet_ldlt(H);

  auto factor_t0 =
      std::chrono::steady_clock::now();

  backend->factorize(H);

  auto factor_t1 =
      std::chrono::steady_clock::now();

  const double factor_ms =
      std::chrono::duration<double,std::milli>(
          factor_t1 - factor_t0).count();

  if (!backend->is_spd()) {
    throw std::runtime_error(
        "Laplace backend reported non-SPD Hessian");
  }

  out.logdet = backend->logdet();

  std::cout << "[factor_ms] "
            << factor_ms
            << "\\n";

  std::cout << "[logdet compare] old="
            << old_logdet
            << " new="
            << out.logdet
            << " diff="
            << (out.logdet - old_logdet)
            << "\\n";"""

if needle not in txt:
    raise RuntimeError("Target block not found.")

txt = txt.replace(needle, replacement)
p.write_text(txt)
print("Patched.")
PY

echo
echo "Diagnostics installed."
