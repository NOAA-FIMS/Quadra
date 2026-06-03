#!/usr/bin/env bash
set -euo pipefail

target="core/laplace/structure_detector.hpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/structure_detector.hpp.exact_precedence.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

small_block = '''    if (options_.prefer_dense_for_small_matrices &&
        info.rows <= options_.dense_size_cutoff) {
      rec.backend = LaplaceBackendKind::DenseLDLT;
      rec.supports_symbolic_reuse = false;
      rec.reason = "small matrix; dense LDLT preferred";
      return rec;
    }

'''

if small_block not in s:
    raise SystemExit("Could not find small-matrix dense block")

s = s.replace(small_block, "", 1)

anchor = '''    if (info.max_bandwidth <= options_.banded_width_cutoff &&
        info.fill_ratio < options_.dense_fill_ratio) {
      rec.backend = LaplaceBackendKind::Banded;
      rec.structure = HessianStructure::Banded;
      rec.supports_symbolic_reuse = true;
      rec.reason = "fixed narrow band detected";
      return rec;
    }

'''

if anchor not in s:
    raise SystemExit("Could not find banded branch anchor")

s = s.replace(anchor, anchor + small_block, 1)
p.write_text(s)
PY

cat <<'EOF'

Reordered StructureDetector rules so exact structure wins before small-dense heuristic.

Run:
  ./run_structure_detector_registry_test.sh
  ./run_laplace_backend_factory_test.sh

Then rebuild and rerun:
  c++ -std=c++17 -O2 -DNDEBUG -g     -Iexternal/Eigen     -Iexternal/LBFGSpp/include     examples/state_space_surplus_production/laplace_state_space_surplus_latent_tridiagonal.cpp     -o build/examples/laplace_state_space_surplus_latent_tridiagonal

  ./build/examples/laplace_state_space_surplus_latent_tridiagonal 1

Expected:
  [recommendation] backend=tridiagonal ... bandwidth=1 ... reason=unit bandwidth

EOF
