#!/usr/bin/env bash
set -euo pipefail

CPP="examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp"
GLOB="examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock_adgraph_global.cpp"

for f in "$CPP" "$GLOB"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: $f not found. Run from repo root."
    exit 1
  fi
done

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$CPP" "${CPP}.before_api_fix.${STAMP}"
cp "$GLOB" "${GLOB}.before_api_fix.${STAMP}"

python3 - <<'PY'
from pathlib import Path
import re

cpp = Path("examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp")
s = cpp.read_text()

old = '''quadra::ParameterVector make_params(std::size_t n_years) {
  quadra::ParameterVector p;
  p.add("log_r0", 8.0, false);
  p.add("log_fbar", -3.7, false);
  p.add("log_q", -7.0, false);
  p.add("logit_sel_a50", -0.2, false);
  p.add("log_sel_slope", 0.0, false);
  for (std::size_t i = 0; i < n_years; ++i)
    p.add("log_rec_dev_" + std::to_string(i + 1), 0.0, true);
  return p;
}
'''

new = '''quadra::ParameterVector make_params(std::size_t n_years) {
  quadra::ParameterVector p;

  auto add_param = [&](const std::string &name, double value, bool random) {
    quadra::Parameter par;
    par.name = name;
    par.value = value;
    par.is_random = random;
    p.add(par);
  };

  add_param("log_r0", 8.0, false);
  add_param("log_fbar", -3.7, false);
  add_param("log_q", -7.0, false);
  add_param("logit_sel_a50", -0.2, false);
  add_param("log_sel_slope", 0.0, false);

  for (std::size_t i = 0; i < n_years; ++i) {
    add_param("log_rec_dev_" + std::to_string(i + 1), 0.0, true);
  }

  return p;
}
'''

if old not in s:
    raise SystemExit("Could not find make_params block to replace.")
s = s.replace(old, new, 1)
cpp.write_text(s)

glob = Path("examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock_adgraph_global.cpp")
g = glob.read_text()
g = g.replace("namespace had { ADGraph *g_ADGraph = nullptr; }",
              "namespace had { threadDefine ADGraph *g_ADGraph = nullptr; }")
glob.write_text(g)

print("Fixed ParameterVector construction and thread-local AD graph global.")
PY

echo
echo "Run:"
echo "  examples/NMFS/afsc_walleye_pollock/run_walleye_pollock_example.sh"
