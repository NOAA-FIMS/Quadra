#!/usr/bin/env bash
set -euo pipefail

L9="examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality"
REPORT="$L9/reports/bigeye_fit_reports.hpp"
SUITE="$L9/reports/bigeye_report_suite.hpp"

if [[ ! -f "$REPORT" || ! -f "$SUITE" ]]; then
  echo "ERROR: missing Level 9 report files. Run from repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$REPORT" "${REPORT}.before_report_compile_fix.${STAMP}"
cp "$SUITE" "${SUITE}.before_report_compile_fix.${STAMP}"

# Fix OptResult member name.
perl -pi -e 's/fit\.joint_value/fit.joint/g' "$REPORT"

# Replace the stale report-suite wrapper with a thin compatibility wrapper.
cat > "$SUITE" <<'EOF_SUITE'
#pragma once

#include "bigeye_fit_reports.hpp"

#include "../objective/bigeye_quadra_objective.hpp"
#include "../../../../../core/optimizer.hpp"

#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

// Compatibility alias for older Level drivers that expect BigeyeReportPaths.
using BigeyeReportPaths = BigeyeFitReportPaths;

inline BigeyeReportPaths default_bigeye_report_paths()
{
  return BigeyeReportPaths{};
}

// Compatibility wrapper for the Level 9 report suite. The observations vector
// is retained in the signature so older driver code does not need to change;
// Level 9 reports use the objective because this level has fleet observations
// and a different fixed-effect parameter ordering.
inline void write_bigeye_report_suite(
    const BigeyeReportPaths &paths,
    const std::vector<Observation> &observations,
    const BigeyeQuadraObjective &objective,
    const quadra::ParameterVector &params,
    const quadra::OptResult &fit)
{
  (void)observations;
  (void)params;

  write_bigeye_fit_summary(paths.summary, fit);
  write_bigeye_fitted_trajectory(paths.trajectory, objective, fit);
  write_bigeye_residual_diagnostics(paths.residual_diagnostics, objective, fit);
  write_bigeye_selectivity_at_age(paths.selectivity, fit);
  write_bigeye_recruitment_deviations(paths.recruitment_deviations, objective,
                                      fit);
  write_bigeye_objective_components(paths.objective_components, objective, fit);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::BigeyeReportPaths;
using pifsc_bigeye_tuna::default_bigeye_report_paths;
using pifsc_bigeye_tuna::write_bigeye_report_suite;
EOF_SUITE

cat > inspect_bigeye_level9_report_compile_fix.sh <<'EOF_INSPECT'
#!/usr/bin/env bash
set -euo pipefail

L9="examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality"
REPORT="$L9/reports/bigeye_fit_reports.hpp"
SUITE="$L9/reports/bigeye_report_suite.hpp"

echo "== OptResult joint member usage =="
grep -n "fit\\.joint\\|fit\\.joint_value" "$REPORT" || true

echo
echo "== Report suite wrapper =="
sed -n '1,120p' "$SUITE"

echo
echo "== Report functions available =="
grep -n "write_bigeye_fit_summary\\|write_bigeye_objective_components\\|write_bigeye_report_suite" "$REPORT" "$SUITE"
EOF_INSPECT

chmod +x inspect_bigeye_level9_report_compile_fix.sh

cat > run_bigeye_level9_report_compile_fix_check.sh <<'EOF_RUN'
#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level9_report_compile_fix.sh

echo
echo "== O3 build Bigeye Level 9 report compile fix =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/quadra/bigeye_level9_estimated_natural_mortality.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level9_report_compile_fix_check

echo
echo "== Run Bigeye Level 9 report compile fix =="
./build/examples/pifsc_bigeye_level9_report_compile_fix_check

echo
echo "== Objective consistency check =="
cat examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_consistency_check.csv

echo
echo "== Patched objective components =="
cat examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_components.csv

echo
echo "== Consistency comparison =="
python3 - <<'PY'
import csv
from pathlib import Path

cons = Path("examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_consistency_check.csv")
comp = Path("examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_components.csv")

vals = {}
with cons.open() as f:
    for r in csv.DictReader(f):
        vals[r["metric"]] = float(r["value"])

components = {}
with comp.open() as f:
    for r in csv.DictReader(f):
        components[r["component"]] = float(r["value"])

print(f'direct_joint_objective={vals["direct_joint_objective"]:.12f}')
print(f'component_joint_total={components["joint_total"]:.12f}')
print(f'difference={components["joint_total"] - vals["direct_joint_objective"]:.12g}')
PY
EOF_RUN

chmod +x run_bigeye_level9_report_compile_fix_check.sh

echo "Patched Level 9 report compile errors."
echo
echo "Run:"
echo "  ./inspect_bigeye_level9_report_compile_fix.sh"
echo "  ./run_bigeye_level9_report_compile_fix_check.sh"
