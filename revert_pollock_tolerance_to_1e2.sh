#!/usr/bin/env bash
set -euo pipefail

for fname in run_pollock_showcase_report.sh run_pollock_polished_report.sh; do
  [[ -f "$fname" ]] || continue
  python3 - "$fname" <<'PYCODE'
from pathlib import Path
import sys
p = Path(sys.argv[1])
s = p.read_text()
s = s.replace("-DQUADRA_LBFGS_GRAD_TOL=1.0e-4", "-DQUADRA_LBFGS_GRAD_TOL=1.0e-2")
s = s.replace("QUADRA_LBFGS_GRAD_TOL=1.0e-4", "QUADRA_LBFGS_GRAD_TOL=1.0e-2")
p.write_text(s)
PYCODE
done

echo "Reverted Pollock showcase/polished runners to QUADRA_LBFGS_GRAD_TOL=1.0e-2."
