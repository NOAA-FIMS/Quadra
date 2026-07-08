#!/usr/bin/env bash
set -euo pipefail

L17="examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic"

if [[ ! -d "$L17" ]]; then
  echo "ERROR: missing $L17."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
for f in \
  "$L17/objective/bigeye_quadra_objective.hpp" \
  "$L17/reports/bigeye_fit_reports.hpp" \
  "$L17/diagnostics/bigeye_initial_numbers_diagnostics.hpp" \
  "$L17/diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
do
  [[ -f "$f" ]] && cp "$f" "${f}.before_layout_fix.${STAMP}"
done

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic")

# Actual driver layout from inspection:
# 0 log_r0
# 1 log_fbar
# 2 log_q_purse_seine
# 3 logit_sel_a50_longline
# 4 log_sel_slope_longline
# 5..14 init_log_number_dev_age_1..10
# 15 log_juvenile_m_multiplier
# 16..25 logit_sel_purse_seine_age_1..10
# 26.. recruitment deviations

# Objective: align layout and move sigma before use.
p = root / "objective/bigeye_quadra_objective.hpp"
s = p.read_text()

s = s.replace("constexpr int kBaseFixed = 6;", "constexpr int kBaseFixed = 5;")
s = s.replace(
    "constexpr int kPurseSeineSelOffset = kBaseFixed + kInitialDevs;",
    "constexpr int kJuvenileMIndex = kBaseFixed + kInitialDevs;\n    constexpr int kPurseSeineSelOffset = kJuvenileMIndex + 1;"
)
s = s.replace(
    "const T log_juvenile_m_multiplier = par[5];",
    "const T log_juvenile_m_multiplier = par[kJuvenileMIndex];"
)

# Remove any duplicate late sigma declaration, then add it right with the other sigmas.
s = re.sub(r"\n\s*const T sigma_log_juvenile_m_multiplier = T\(0\.50\);", "", s)
s = s.replace(
    "const T sigma_ps_sel_dev = T(1.0);",
    "const T sigma_ps_sel_dev = T(1.0);\n    const T sigma_log_juvenile_m_multiplier = T(0.50);",
    1
)

# Ensure juvenile prior is after sigma declaration and after base priors.
prior = "    nll = nll + T(0.5) * square_t(log_juvenile_m_multiplier / sigma_log_juvenile_m_multiplier);\n"
s = s.replace(prior, "")
anchor = "    nll = nll + normal_prior(log_sel_slope_longline, std::log(1.2), 0.35);\n"
if anchor not in s:
    raise SystemExit("Could not find slope prior anchor in objective")
s = s.replace(anchor, anchor + prior, 1)

p.write_text(s)

# Reports: initial numbers at 5+a; juvenile at 15; PS selectivity at 16+a.
p = root / "reports/bigeye_fit_reports.hpp"
if p.exists():
    s = p.read_text()

    s = s.replace("const int idx = 6 + a;", "const int idx = 5 + a;")
    s = s.replace(
        "const int idx = 6 + pifsc_bigeye_tuna::kAges + a;",
        "const int idx = 5 + pifsc_bigeye_tuna::kAges + 1 + a;"
    )

    # Replace juvenile report indices from fit.par[5] to layout index.
    s = s.replace('fit.par[5] << "\\n";', 'fit.par[5 + pifsc_bigeye_tuna::kAges] << "\\n";')
    s = s.replace('std::exp(fit.par[5]) << "\\n";', 'std::exp(fit.par[5 + pifsc_bigeye_tuna::kAges]) << "\\n";')
    s = s.replace('0.45 * std::exp(fit.par[5]) << "\\n";', '0.45 * std::exp(fit.par[5 + pifsc_bigeye_tuna::kAges]) << "\\n";')

    p.write_text(s)

# Initial numbers diagnostic: baseFixed=5, juvenile index=15, no duplicate stale M.
p = root / "diagnostics/bigeye_initial_numbers_diagnostics.hpp"
if p.exists():
    s = p.read_text()
    s = s.replace("constexpr int kBaseFixed = 6;", "constexpr int kBaseFixed = 5;")
    s = s.replace("std::exp(fit.par[5])", "std::exp(fit.par[kBaseFixed + kAges])")
    s = re.sub(r"\n\s*const double log_m = std::log\(0\.45\);", "", s)
    s = re.sub(r"\n\s*const double m = std::exp\(log_m\);", "", s)
    p.write_text(s)

# Age-comp residual diagnostic: baseFixed=5, juvenile index=15, PS offset=16.
p = root / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
if p.exists():
    s = p.read_text()
    s = s.replace("constexpr int kBaseFixed = 6;", "constexpr int kBaseFixed = 5;")
    s = s.replace(
        "constexpr int kPurseSeineSelOffset = kBaseFixed + kInitialDevs;",
        "constexpr int kJuvenileMIndex = kBaseFixed + kInitialDevs;\n  constexpr int kPurseSeineSelOffset = kJuvenileMIndex + 1;"
    )
    s = s.replace("std::exp(fit.par[5])", "std::exp(fit.par[kJuvenileMIndex])")
    p.write_text(s)

PY

echo "Aligned Level 17 parameter layout to actual driver order."
echo
echo "Expected layout:"
echo "  0..4   base fixed effects"
echo "  5..14  initial number deviations"
echo "  15     log_juvenile_m_multiplier"
echo "  16..25 purse-seine age selectivity logits"
echo "  26..   recruitment deviations"
echo
echo "Run:"
echo "  ./inspect_bigeye_level17_juvenile_mortality.sh"
echo "  ./run_bigeye_level17_juvenile_mortality_check.sh"
