#!/usr/bin/env bash
set -euo pipefail

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_state_space_latent_runtime_phase2_v2.sh first."
  exit 1
fi

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()
s = s.replace("ss::make_example_data()", "ss::make_demo_data()")
s = s.replace("ss::make_example_parameters()", "ss::make_demo_parameters()")
p.write_text(s)
PY

echo "Repaired Phase 2 runtime source: make_example_* -> make_demo_*"
echo
echo "Run:"
echo "  ./run_state_space_surplus_latent_runtime_phase2.sh 20"
