#!/usr/bin/env bash
set -euo pipefail

# repair_state_space_latent_runtime_derivative_cache_output_v1.sh
#
# Repairs partial derivative-cache install by adding only the missing
# derivative_cached timing and output lines.
#
# Assumes install_state_space_latent_runtime_derivative_cache_v1.sh already
# added evaluate_cached_derivatives() and derivative cache members.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.derivative_cache_output.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if "evaluate_cached_derivatives()" not in s:
    raise SystemExit("evaluate_cached_derivatives() not found; rerun derivative-cache installer first")

if "derivative_cached_total_ms" not in s:
    anchor = "  const auto result_cached0 = Clock::now();\n"
    insert = (
        "  const auto derivative_cached0 = Clock::now();\n"
        "  for (int r = 0; r < reps; ++r) {\n"
        "    last = runtime.evaluate_cached_derivatives();\n"
        "  }\n"
        "  const auto derivative_cached1 = Clock::now();\n"
        "  const double derivative_cached_total_ms =\n"
        "      ms_between(derivative_cached0, derivative_cached1);\n"
        "  const double derivative_cached_avg_ms =\n"
        "      derivative_cached_total_ms / static_cast<double>(reps);\n\n"
    )
    if anchor not in s:
        raise SystemExit("Could not find result_cached0 timing anchor")
    s = s.replace(anchor, insert + anchor, 1)

if 'derivative_cached_avg_ms = ' not in s:
    anchor = '  std::cout << "result_cached_total_ms = " << result_cached_total_ms << "\\n";\n'
    insert = (
        '  std::cout << "derivative_cached_total_ms = " << derivative_cached_total_ms << "\\n";\n'
        '  std::cout << "derivative_cached_avg_ms = " << derivative_cached_avg_ms << "\\n";\n'
    )
    if anchor not in s:
        # Fallback anchor: print before return 0
        anchor2 = "\n  return 0;"
        insert2 = (
            '  std::cout << "derivative_cached_total_ms = " << derivative_cached_total_ms << "\\n";\n'
            '  std::cout << "derivative_cached_avg_ms = " << derivative_cached_avg_ms << "\\n";\n'
        )
        if anchor2 not in s:
            raise SystemExit("Could not find output anchor")
        s = s.replace(anchor2, "\n" + insert2 + anchor2, 1)
    else:
        s = s.replace(anchor, insert + anchor, 1)

p.write_text(s)
print("Installed derivative-cache timing/output repair.")
PY

cat <<'EOF'

Installed derivative-cache timing/output repair.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  derivative_cached_total_ms = ...
  derivative_cached_avg_ms = ...

EOF
