#!/usr/bin/env bash
set -euo pipefail

L7="examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity"
DRIVER="$L7/quadra/bigeye_level7_dual_age_selectivity.cpp"

if [[ ! -d "$L7" ]]; then
  echo "ERROR: missing $L7. Run from repo root after Level 7 patch."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$DRIVER" "${DRIVER}.before_minimal_report_fix.${STAMP}"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/quadra/bigeye_level7_dual_age_selectivity.cpp")
s = p.read_text()

patterns = [
    r'\n\s*pifsc_bigeye_tuna::write_bigeye_report_suite\([\s\S]*?\);\n',
    r'\n\s*pifsc_bigeye_tuna::write_bigeye_functional_analysis_report\([\s\S]*?\);\n',
    r'\n\s*pifsc_bigeye_tuna::write_bigeye_fixed_effect_geometry_report\([\s\S]*?\);\n',
    r'\n\s*pifsc_bigeye_tuna::write_longline_slope_geometry_scan\([\s\S]*?\);\n',
]
for pat in patterns:
    s = re.sub(pat, "\n", s, count=1)

if "#include <filesystem>" not in s:
    s = s.replace("#include <fstream>\n", "#include <filesystem>\n#include <fstream>\n", 1)

ensure = '    std::filesystem::create_directories("examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs");\n'
if ensure not in s:
    marker = "    auto fit = quadra::optimize_lbfgs(objective, params, opts);\n"
    if marker not in s:
        raise SystemExit("Could not find optimize marker.")
    s = s.replace(marker, marker + "\n" + ensure, 1)

fit_summary = '''
    {
      std::ofstream out("examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_fit_summary.csv");
      if (!out) {
        throw std::runtime_error("Could not open Level 7 fit summary CSV");
      }
      out << std::setprecision(12);
      out << "field,value\\n";
      out << "objective," << fit.value << "\\n";
      out << "grad_norm," << fit.grad_norm << "\\n";
      out << "iterations," << fit.iterations << "\\n";
      out << "converged," << (fit.converged ? "yes" : "no") << "\\n";
      out << "message," << fit.message << "\\n";
      out << "laplace,yes\\n";
      out << "random_effects," << fit.u_hat.size() << "\\n";
      if (fit.par.size() >= 4) {
        out << "log_r0," << fit.par[0] << "\\n";
        out << "r0," << std::exp(fit.par[0]) << "\\n";
        out << "log_fbar," << fit.par[1] << "\\n";
        out << "fbar," << std::exp(fit.par[1]) << "\\n";
        out << "log_q_longline," << fit.par[2] << "\\n";
        out << "q_longline," << std::exp(fit.par[2]) << "\\n";
        out << "log_q_purse_seine," << fit.par[3] << "\\n";
        out << "q_purse_seine," << std::exp(fit.par[3]) << "\\n";
      }
    }

'''

if "bigeye_level7_fit_summary.csv" not in s:
    s = s.replace(ensure, ensure + fit_summary, 1)

print_block = '''
    std::cout << "PIFSC bigeye-tuna-style Level 7 dual-age-selectivity Quadra Laplace recruitment-deviation fit\\n";
    std::cout << "objective:  " << fit.value << "\\n";
    std::cout << "grad_norm:  " << fit.grad_norm << "\\n";
    std::cout << "converged:  " << (fit.converged ? "yes" : "no") << "\\n";
    std::cout << "message:    " << fit.message << "\\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_fit_summary.csv\\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_recruitment_diagnostics.txt\\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_recruitment_diagnostics.csv\\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_safe_fixed_effect_wiggle_diagnostics.txt\\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_safe_fixed_effect_wiggle_diagnostics.csv\\n";
'''
s = re.sub(
    r'\n\s*std::cout\s*<<\s*"PIFSC bigeye[\s\S]*?\n\s*return 0;',
    "\n" + print_block + "\n    return 0;",
    s,
    count=1,
)

p.write_text(s)
PY

cat > inspect_bigeye_level7_minimal_report_fix.sh <<'EOF_INSPECT'
#!/usr/bin/env bash
set -euo pipefail

L7="examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity"
DRIVER="$L7/quadra/bigeye_level7_dual_age_selectivity.cpp"

echo "== Remaining incompatible report calls? =="
grep -n "write_bigeye_report_suite\\|write_bigeye_functional_analysis_report\\|write_bigeye_fixed_effect_geometry_report\\|write_longline_slope_geometry_scan" "$DRIVER" || true

echo
echo "== Level 7 retained diagnostic calls =="
grep -n "write_recruitment_diagnostics\\|write_safe_fixed_effect_wiggle" "$DRIVER" || true

echo
echo "== Level 7 fixed effects =="
grep -n "params.add" "$DRIVER" | head -40
EOF_INSPECT

chmod +x inspect_bigeye_level7_minimal_report_fix.sh

cat > run_bigeye_level7_minimal_report_fix_check.sh <<'EOF_RUN'
#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level7_minimal_report_fix.sh

echo
echo "== O3 build Bigeye Level 7 minimal report fix =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/quadra/bigeye_level7_dual_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level7_minimal_report_fix_check

echo
echo "== Run Bigeye Level 7 minimal report fix =="
./build/examples/pifsc_bigeye_level7_minimal_report_fix_check

echo
echo "== Level 7 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_fit_summary.csv

echo
echo "== Level 7 recruitment diagnostics preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_recruitment_diagnostics.txt

echo
echo "== Level 6 vs Level 7 recruitment comparison =="
python3 - <<'PY'
import csv
from pathlib import Path

base = Path("examples/NMFS/pifsc_bigeye_tuna")
paths = {
    "level6": base / "level6_purse_seine_age_selectivity/outputs/bigeye_level6_recruitment_diagnostics.csv",
    "level7": base / "level7_dual_age_selectivity/outputs/bigeye_level7_recruitment_diagnostics.csv",
}
metrics = ["sd", "lag1_correlation", "roughness", "total_prior_nll", "max_abs"]
print("level," + ",".join(metrics))
for level, p in paths.items():
    vals = {}
    if p.exists():
        with p.open() as f:
            for r in csv.DictReader(f):
                if r["section"] == "summary" and r["metric"] in metrics:
                    vals[r["metric"]] = r["value"]
    print(level + "," + ",".join(vals.get(m, "") for m in metrics))
PY
EOF_RUN

chmod +x run_bigeye_level7_minimal_report_fix_check.sh

echo "Patched Level 7 to bypass incompatible copied report writers and retain recruitment/safe wiggle diagnostics."
echo
echo "Run:"
echo "  ./inspect_bigeye_level7_minimal_report_fix.sh"
echo "  ./run_bigeye_level7_minimal_report_fix_check.sh"
