#!/usr/bin/env bash
set -euo pipefail

LAPLACE="core/laplace.hpp"
DIAG="core/laplace/laplace_gradient_diagnostics.hpp"

if [[ ! -f "$LAPLACE" ]]; then
  echo "ERROR: $LAPLACE not found. Run this from the Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${LAPLACE}.before_diagnostics_header_cleanup.${STAMP}"
cp "$LAPLACE" "$BACKUP"
echo "Backed up $LAPLACE to $BACKUP"

mkdir -p "$(dirname "$DIAG")"

python3 - <<'PY'
from pathlib import Path
import re

diag = Path("core/laplace/laplace_gradient_diagnostics.hpp")
diag.write_text(r'''#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <iostream>

namespace quadra {
namespace laplace {
namespace diagnostics {

inline void print_du_dtheta_summary(const Eigen::MatrixXd &dU) {
#ifdef QUADRA_DEBUG_DU_DTHETA_NORMS
  std::cout << "Quadra dU diagnostic\n";

  std::cout << "  dU_col_norms = ";
  for (Eigen::Index j = 0; j < dU.cols(); ++j) {
    std::cout << dU.col(j).norm();
    if (j + 1 < dU.cols()) std::cout << " ";
  }
  std::cout << "\n";

  std::cout << "  dU_col_maxabs = ";
  for (Eigen::Index j = 0; j < dU.cols(); ++j) {
    std::cout << dU.col(j).cwiseAbs().maxCoeff();
    if (j + 1 < dU.cols()) std::cout << " ";
  }
  std::cout << "\n";

  std::cout << "  dU_first_rows =";
  const Eigen::Index nprint = std::min<Eigen::Index>(5, dU.rows());
  for (Eigen::Index r = 0; r < nprint; ++r) {
    std::cout << "\n    row " << r << ": ";
    for (Eigen::Index j = 0; j < dU.cols(); ++j) {
      std::cout << dU(r, j);
      if (j + 1 < dU.cols()) std::cout << " ";
    }
  }
  std::cout << "\n";
#else
  (void)dU;
#endif
}

inline void print_theta_only_vs_total_logdet_gradient(
    const Eigen::VectorXd &theta_only, const Eigen::VectorXd &total) {
#ifdef QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL
  std::cout << "Quadra logdet Hdot diagnostic\n";
  std::cout << "  theta_only_logdet_grad = " << theta_only.transpose()
            << "\n";
  std::cout << "  total_logdet_grad      = " << total.transpose() << "\n";
  std::cout << "  implicit_u_contribution= "
            << (total - theta_only).transpose() << "\n";
#else
  (void)theta_only;
  (void)total;
#endif
}

inline void print_hdot_exact_vs_fd_trace(
    const Eigen::VectorXd &exact_trace, const Eigen::VectorXd &fd_trace,
    const Eigen::VectorXd &rel_hdot_matrix_err) {
#ifdef QUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE
  std::cout << "Quadra Hdot exact-vs-FD trace diagnostic\n";
  std::cout << "  exact_total_logdet_grad = "
            << exact_trace.transpose() << "\n";
  std::cout << "  fd_total_logdet_grad    = " << fd_trace.transpose()
            << "\n";
  std::cout << "  exact_minus_fd          = "
            << (exact_trace - fd_trace).transpose() << "\n";
  std::cout << "  rel_Hdot_matrix_err     = "
            << rel_hdot_matrix_err.transpose() << "\n";
#else
  (void)exact_trace;
  (void)fd_trace;
  (void)rel_hdot_matrix_err;
#endif
}

inline void print_gradient_parts(const Eigen::VectorXd &joint_grad,
                                 const Eigen::VectorXd &logdet_grad,
                                 const Eigen::VectorXd &total_grad) {
#ifdef QUADRA_DEBUG_LAPLACE_GRADIENT_PARTS
  std::cout << "Quadra gradient parts\n";
  std::cout << "  joint_grad  = " << joint_grad.transpose() << "\n";
  std::cout << "  logdet_grad = " << logdet_grad.transpose() << "\n";
  std::cout << "  total_grad  = " << total_grad.transpose() << "\n";
#else
  (void)joint_grad;
  (void)logdet_grad;
  (void)total_grad;
#endif
}

inline void print_logdet_gradient_comparison(
    const Eigen::VectorXd &exact_logdet_grad,
    const Eigen::VectorXd &fd_logdet_grad) {
#ifdef QUADRA_DEBUG_LAPLACE_GRADIENT_PARTS
  std::cout << "Quadra logdet gradient parts\n";
  std::cout << "  logdet_grad      = " << exact_logdet_grad.transpose()
            << "\n";
  std::cout << "  logdet_fd_grad   = " << fd_logdet_grad.transpose()
            << "\n";
  std::cout << "  logdet_grad diff = "
            << (exact_logdet_grad - fd_logdet_grad).transpose() << "\n";
#else
  (void)exact_logdet_grad;
  (void)fd_logdet_grad;
#endif
}

}  // namespace diagnostics
}  // namespace laplace
}  // namespace quadra
''')

path = Path("core/laplace.hpp")
text = path.read_text()

include_line = '#include "laplace/laplace_gradient_diagnostics.hpp"\n'
if include_line not in text:
    marker = '#include "laplace/had_quadra_replay_reuse_sparse_hdot_provider.hpp"\n'
    if marker in text:
        text = text.replace(marker, marker + include_line, 1)
    else:
        matches = list(re.finditer(r'^#include .*$\n?', text, re.M))
        if matches:
            pos = matches[-1].end()
            text = text[:pos] + include_line + text[pos:]
        else:
            text = include_line + text

# Replace temporary dU inline diagnostic with helper call.
du_block = re.compile(
    r'\n#ifdef QUADRA_DEBUG_DU_DTHETA_NORMS\n'
    r'  \{\n'
    r'    std::cout << "Quadra dU diagnostic\\n";.*?'
    r'  \}\n'
    r'#endif\n',
    re.S,
)
text, n_du = du_block.subn(
    '\n  laplace::diagnostics::print_du_dtheta_summary(dU);\n',
    text,
)

du_assign = (
    '  Eigen::MatrixXd dU =\n'
    '      implicit_du_dtheta_all(model, params, theta, u_hat, &H_factor, &solver);\n\n'
)
if n_du == 0 and du_assign in text and "print_du_dtheta_summary(dU)" not in text:
    text = text.replace(
        du_assign,
        du_assign + '  laplace::diagnostics::print_du_dtheta_summary(dU);\n\n',
        1,
    )

# Replace theta-only print section if a temporary inline block exists.
theta_block = re.compile(
    r'\n#ifdef QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL\n'
    r'  \{\n'
    r'    const Eigen::MatrixXd zero_dU =.*?'
    r'    std::cout << "Quadra logdet Hdot diagnostic\\n";.*?'
    r'  \}\n'
    r'#endif\n',
    re.S,
)
theta_replacement = (
    '\n#ifdef QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL\n'
    '  {\n'
    '    const Eigen::MatrixXd zero_dU =\n'
    '        Eigen::MatrixXd::Zero(u_hat.size(), theta.size());\n\n'
    '    const auto Hdots_theta_only = random_hessian_directional_exact_all(\n'
    '        model, params, theta, u_hat, zero_dU, get_pattern_for_logdet);\n\n'
    '    Eigen::VectorXd theta_only = Eigen::VectorXd::Zero(theta.size());\n'
    '    for (Eigen::Index i = 0; i < theta.size(); ++i) {\n'
    '      theta_only[i] =\n'
    '          0.5 * logdet_directional_derivative_from_hdot(\n'
    '                    solver, Hdots_theta_only[static_cast<std::size_t>(i)],\n'
    '                    options);\n'
    '    }\n\n'
    '    laplace::diagnostics::print_theta_only_vs_total_logdet_gradient(\n'
    '        theta_only, grad);\n'
    '  }\n'
    '#endif\n'
)
text, _ = theta_block.subn(theta_replacement, text)

# Replace exact-vs-FD print section if a temporary inline block exists.
hdot_block = re.compile(
    r'\n#ifdef QUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE\n'
    r'  \{\n'
    r'    Eigen::VectorXd fd_trace =.*?'
    r'    std::cout << "Quadra Hdot exact-vs-FD trace diagnostic\\n";.*?'
    r'  \}\n'
    r'#endif\n',
    re.S,
)
hdot_replacement = (
    '\n#ifdef QUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE\n'
    '  {\n'
    '    Eigen::VectorXd fd_trace = Eigen::VectorXd::Zero(theta.size());\n'
    '    Eigen::VectorXd exact_trace = Eigen::VectorXd::Zero(theta.size());\n'
    '    Eigen::VectorXd rel_hdot_err = Eigen::VectorXd::Zero(theta.size());\n\n'
    '    for (Eigen::Index i = 0; i < theta.size(); ++i) {\n'
    '      const Eigen::SparseMatrix<double> Hdot_fd =\n'
    '          random_hessian_directional_implicit_fd_with_du(\n'
    '              model, params, theta, u_hat, i, dU.col(i), 1.0e-5);\n\n'
    '      const Eigen::SparseMatrix<double> &Hdot_exact =\n'
    '          Hdots[static_cast<std::size_t>(i)];\n\n'
    '      fd_trace[i] =\n'
    '          0.5 * logdet_directional_derivative_from_hdot(\n'
    '                    solver, Hdot_fd, options);\n'
    '      exact_trace[i] =\n'
    '          0.5 * logdet_directional_derivative_from_hdot(\n'
    '                    solver, Hdot_exact, options);\n\n'
    '      const Eigen::SparseMatrix<double> diff = Hdot_exact - Hdot_fd;\n'
    '      rel_hdot_err[i] =\n'
    '          diff.norm() / std::max(1.0e-12, Hdot_fd.norm());\n'
    '    }\n\n'
    '    laplace::diagnostics::print_hdot_exact_vs_fd_trace(\n'
    '        exact_trace, fd_trace, rel_hdot_err);\n'
    '  }\n'
    '#endif\n'
)
text, _ = hdot_block.subn(hdot_replacement, text)

path.write_text(text)
PY

echo
echo "Created diagnostics header:"
echo "  $DIAG"
echo
echo "Remaining diagnostic macro references:"
grep -R "QUADRA_DEBUG_LOGDET_THETA_ONLY_VS_TOTAL\|QUADRA_DEBUG_HDOT_EXACT_VS_FD_TRACE\|QUADRA_DEBUG_DU_DTHETA_NORMS\|QUADRA_DEBUG_LAPLACE_GRADIENT_PARTS" core/laplace.hpp "$DIAG" || true
echo
echo "Suggested clean build:"
echo 'clang++ -std=c++17 -g -I"external/eigen/" examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp'
