#!/usr/bin/env bash
set -euo pipefail

FILE="core/laplace.hpp"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run from Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${FILE}.before_logdet_theta_only_diagnostic.${STAMP}"
cp "$FILE" "$BACKUP"
echo "Backed up $FILE to:"
echo "  $BACKUP"

python3 - <<'PY'
from pathlib import Path

path = Path("core/laplace.hpp")
text = path.read_text()

if "QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL" in text:
    print("Diagnostic already installed.")
    raise SystemExit(0)

needle = """  const auto Hdots = random_hessian_directional_exact_all(
      model, params, theta, u_hat, dU, get_pattern_for_logdet);

  for (Eigen::Index i = 0; i < theta.size(); ++i) {
    const Eigen::SparseMatrix<double> &Hdot =
        Hdots[static_cast<std::size_t>(i)];

    grad[i] =
        0.5 * logdet_directional_derivative_from_hdot(solver, Hdot, options);
  }
"""

replacement = """  const auto Hdots = random_hessian_directional_exact_all(
      model, params, theta, u_hat, dU, get_pattern_for_logdet);

  for (Eigen::Index i = 0; i < theta.size(); ++i) {
    const Eigen::SparseMatrix<double> &Hdot =
        Hdots[static_cast<std::size_t>(i)];

    grad[i] =
        0.5 * logdet_directional_derivative_from_hdot(solver, Hdot, options);
  }

#ifdef QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL
  {
    const Eigen::MatrixXd zero_dU =
        Eigen::MatrixXd::Zero(u_hat.size(), theta.size());

    const auto Hdots_theta_only = random_hessian_directional_exact_all(
        model, params, theta, u_hat, zero_dU, get_pattern_for_logdet);

    Eigen::VectorXd theta_only = Eigen::VectorXd::Zero(theta.size());
    for (Eigen::Index i = 0; i < theta.size(); ++i) {
      theta_only[i] =
          0.5 * logdet_directional_derivative_from_hdot(
                    solver, Hdots_theta_only[static_cast<std::size_t>(i)],
                    options);
    }

    std::cout << "Quadra logdet Hdot diagnostic\\n";
    std::cout << "  theta_only_logdet_grad = "
              << theta_only.transpose() << "\\n";
    std::cout << "  total_logdet_grad      = "
              << grad.transpose() << "\\n";
    std::cout << "  implicit_u_contribution= "
              << (grad - theta_only).transpose() << "\\n";
  }
#endif
"""

if needle not in text:
    raise RuntimeError("Could not find target Hdot loop in core/laplace.hpp")

text = text.replace(needle, replacement, 1)
path.write_text(text)
print("Installed QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL diagnostic.")
PY

echo
echo "Build with:"
echo 'clang++ -std=c++17 -g -I"external/eigen/" -DQUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp'
