#!/usr/bin/env bash
set -euo pipefail

FILE="core/laplace.hpp"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run from Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${FILE}.before_hdot_exact_vs_fd_trace_diagnostic.${STAMP}"
cp "$FILE" "$BACKUP"
echo "Backed up $FILE to:"
echo "  $BACKUP"

python3 - <<'PY'
from pathlib import Path

path = Path("core/laplace.hpp")
text = path.read_text()

if "QUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE" in text:
    print("Diagnostic already installed.")
    raise SystemExit(0)

needle = """#ifdef QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL
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

insert = needle + """

#ifdef QUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE
  {
    Eigen::VectorXd fd_trace = Eigen::VectorXd::Zero(theta.size());
    Eigen::VectorXd exact_trace = Eigen::VectorXd::Zero(theta.size());
    Eigen::VectorXd rel_hdot_err = Eigen::VectorXd::Zero(theta.size());

    for (Eigen::Index i = 0; i < theta.size(); ++i) {
      const Eigen::SparseMatrix<double> Hdot_fd =
          random_hessian_directional_implicit_fd_with_du(
              model, params, theta, u_hat, i, dU.col(i), 1.0e-5);

      const Eigen::SparseMatrix<double> &Hdot_exact =
          Hdots[static_cast<std::size_t>(i)];

      fd_trace[i] =
          0.5 * logdet_directional_derivative_from_hdot(
                    solver, Hdot_fd, options);
      exact_trace[i] =
          0.5 * logdet_directional_derivative_from_hdot(
                    solver, Hdot_exact, options);

      const Eigen::SparseMatrix<double> diff = Hdot_exact - Hdot_fd;
      rel_hdot_err[i] =
          diff.norm() / std::max(1.0e-12, Hdot_fd.norm());
    }

    std::cout << "Quadra Hdot exact-vs-FD trace diagnostic\\n";
    std::cout << "  exact_total_logdet_grad = "
              << exact_trace.transpose() << "\\n";
    std::cout << "  fd_total_logdet_grad    = "
              << fd_trace.transpose() << "\\n";
    std::cout << "  exact_minus_fd          = "
              << (exact_trace - fd_trace).transpose() << "\\n";
    std::cout << "  rel_Hdot_matrix_err     = "
              << rel_hdot_err.transpose() << "\\n";
  }
#endif
"""

if needle not in text:
    raise RuntimeError("Could not find theta-only diagnostic block. Install that first or restore target structure.")

text = text.replace(needle, insert, 1)
path.write_text(text)
print("Installed QUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE diagnostic.")
PY

echo
echo "Build with:"
echo 'clang++ -std=c++17 -g -I"external/eigen/" -DQUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp'
