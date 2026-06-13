#!/usr/bin/env bash
set -euo pipefail

FILE="core/laplace.hpp"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run from Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${FILE}.before_du_dtheta_norm_diagnostic.${STAMP}"
cp "$FILE" "$BACKUP"
echo "Backed up $FILE to:"
echo "  $BACKUP"

python3 - <<'PY'
from pathlib import Path

path = Path("core/laplace.hpp")
text = path.read_text()

if "QUADRA_DEBUG_DU_DTHETA_NORMS" in text:
    print("dU norm diagnostic already installed.")
    raise SystemExit(0)

needle = """  Eigen::MatrixXd dU =
      implicit_du_dtheta_all(model, params, theta, u_hat, &H_factor, &solver);

  const auto timing_du_end = std::chrono::steady_clock::now();
"""

replacement = """  Eigen::MatrixXd dU =
      implicit_du_dtheta_all(model, params, theta, u_hat, &H_factor, &solver);

#ifdef QUADRA_DEBUG_DU_DTHETA_NORMS
  {
    std::cout << "Quadra dU diagnostic\\n";
    std::cout << "  dU_col_norms = ";
    for (Eigen::Index j = 0; j < dU.cols(); ++j) {
      std::cout << dU.col(j).norm();
      if (j + 1 < dU.cols()) {
        std::cout << " ";
      }
    }
    std::cout << "\\n";

    std::cout << "  dU_col_maxabs = ";
    for (Eigen::Index j = 0; j < dU.cols(); ++j) {
      std::cout << dU.col(j).cwiseAbs().maxCoeff();
      if (j + 1 < dU.cols()) {
        std::cout << " ";
      }
    }
    std::cout << "\\n";

    std::cout << "  dU_first_rows =";
    const Eigen::Index nprint = std::min<Eigen::Index>(5, dU.rows());
    for (Eigen::Index r = 0; r < nprint; ++r) {
      std::cout << "\\n    row " << r << ": ";
      for (Eigen::Index j = 0; j < dU.cols(); ++j) {
        std::cout << dU(r, j);
        if (j + 1 < dU.cols()) {
          std::cout << " ";
        }
      }
    }
    std::cout << "\\n";
  }
#endif

  const auto timing_du_end = std::chrono::steady_clock::now();
"""

if needle not in text:
    raise RuntimeError("Could not find dU assignment block in core/laplace.hpp")

text = text.replace(needle, replacement, 1)
path.write_text(text)
print("Installed QUADRA_DEBUG_DU_DTHETA_NORMS diagnostic.")
PY

echo
echo "Build with:"
echo 'clang++ -std=c++17 -g -I"external/eigen/" -DQUADRA_DEBUG_DU_DTHETA_NORMS examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp'
