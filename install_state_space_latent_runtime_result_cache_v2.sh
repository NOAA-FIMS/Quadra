#!/usr/bin/env bash
set -euo pipefail

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.result_cache_v2.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if 'evaluate_cached_result()' in s and 'result_cached_avg_ms' in s:
    print('Result cache patch already appears installed.')
    raise SystemExit(0)

if 'bool result_cached_' not in s:
    old = '  bool initialized_;\n\n  quadra::laplace::BackendRecommendation recommendation_;'
    new = '  bool initialized_;\n  bool result_cached_ = false;\n  EvalResult cached_result_;\n\n  quadra::laplace::BackendRecommendation recommendation_;'
    if old not in s:
        raise SystemExit('Could not find runtime member anchor')
    s = s.replace(old, new, 1)

if 'cached_result_ = out;' not in s:
    old = '    out.objective = out.joint + out.correction;\n\n    return out;\n  }'
    new = '    out.objective = out.joint + out.correction;\n\n    cached_result_ = out;\n    result_cached_ = true;\n\n    return out;\n  }'
    if old not in s:
        raise SystemExit('Could not find evaluate_at_xhat return anchor')
    s = s.replace(old, new, 1)

if 'EvalResult evaluate_cached_result() const' not in s:
    anchor = '  const char* backend_name() const {'
    method = (
        '  EvalResult evaluate_cached_result() const {\n'
        '    if (!initialized_ || !result_cached_) {\n'
        '      throw std::runtime_error("runtime result cache used before initialization");\n'
        '    }\n'
        '    return cached_result_;\n'
        '  }\n\n'
    )
    if anchor not in s:
        raise SystemExit('Could not find backend_name anchor')
    s = s.replace(anchor, method + anchor, 1)

if 'result_cached_total_ms' not in s:
    anchor = '  const double cached_avg_ms = cached_total_ms / static_cast<double>(reps);\n'
    insert = (
        '\n  const auto result_cached0 = Clock::now();\n'
        '  for (int r = 0; r < reps; ++r) {\n'
        '    last = runtime.evaluate_cached_result();\n'
        '  }\n'
        '  const auto result_cached1 = Clock::now();\n'
        '  const double result_cached_total_ms =\n'
        '      ms_between(result_cached0, result_cached1);\n'
        '  const double result_cached_avg_ms =\n'
        '      result_cached_total_ms / static_cast<double>(reps);\n'
    )
    if anchor not in s:
        raise SystemExit('Could not find cached_avg_ms timing anchor')
    s = s.replace(anchor, anchor + insert, 1)

if 'result_cached_avg_ms = ' not in s:
    anchor = '  std::cout << "cached_avg_ms = " << cached_avg_ms << "\\n";\n'
    insert = (
        '  std::cout << "result_cached_total_ms = " << result_cached_total_ms << "\\n";\n'
        '  std::cout << "result_cached_avg_ms = " << result_cached_avg_ms << "\\n";\n'
    )
    if anchor not in s:
        raise SystemExit('Could not find cached_avg_ms output anchor')
    s = s.replace(anchor, anchor + insert, 1)

p.write_text(s)
print('Installed robust result-cache patch.')
PY

cat <<'EOF'

Installed robust runtime result-cache benchmark.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected additional output:
  result_cached_total_ms = ...
  result_cached_avg_ms = ...

EOF
