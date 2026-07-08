#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna"

LEVELS=(
  "level20_longline_selectivity_regularization_scan"
  "level21_age_based_natural_mortality_diagnostic"
  "level23_longline_selectivity_smoothness_scan"
)

REPORT="$ROOT/workflow/bigeye_objective_operator_audit.txt"
CSV="$ROOT/workflow/bigeye_objective_operator_audit.csv"

mkdir -p "$ROOT/workflow"

python3 - <<'PY'
from pathlib import Path
import re
import csv

root = Path("examples/NMFS/pifsc_bigeye_tuna")
levels = [
    "level20_longline_selectivity_regularization_scan",
    "level21_age_based_natural_mortality_diagnostic",
    "level23_longline_selectivity_smoothness_scan",
]

report = root / "workflow/bigeye_objective_operator_audit.txt"
csv_path = root / "workflow/bigeye_objective_operator_audit.csv"

checks = []

def add(level, check, status, detail):
    checks.append({
        "level": level,
        "check": check,
        "status": status,
        "detail": detail.replace("\n", " ").strip()
    })

def read(path):
    return path.read_text(errors="replace") if path.exists() else ""

for level in levels:
    base = root / level
    obj = base / "objective/bigeye_quadra_objective.hpp"
    driver = None
    if (base / "quadra").exists():
        drivers = sorted((base / "quadra").glob("bigeye_level*.cpp"))
        driver = drivers[0] if drivers else None
    diag_age = base / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
    diag_init = base / "diagnostics/bigeye_initial_numbers_diagnostics.hpp"
    diag_ll = base / "diagnostics/bigeye_longline_prediction_decomposition.hpp"
    diag_ps = base / "diagnostics/bigeye_purse_seine_prediction_decomposition.hpp"

    s_obj = read(obj)
    s_driver = read(driver) if driver else ""
    s_age = read(diag_age)
    s_init = read(diag_init)
    s_ll = read(diag_ll)
    s_ps = read(diag_ps)

    driver_params = re.findall(r'params\.add\(\{"([^"]+)"', s_driver)
    par_refs = sorted(set(int(x) for x in re.findall(r'par\[(\d+)\]', s_obj)))
    offset_names = sorted(set(re.findall(r'constexpr\s+int\s+(k[A-Za-z0-9_]+)\s*=', s_obj)))
    if driver_params and par_refs:
        add(level, "1_parameter_layout", "inspect",
            f"driver_params={len(driver_params)} direct_par_refs={par_refs} offsets={offset_names}")
    else:
        add(level, "1_parameter_layout", "warn",
            "Could not extract driver params or direct par refs")

    pred_selected = "pred_age_comp[i] = n[i] * sel[i]" in s_obj or "pred_age_comp[i] = n[i] * selectivity[i]" in s_obj
    around_pred = s_obj[max(0, s_obj.find("pred_age_comp")-500):s_obj.find("pred_age_comp")+1200] if "pred_age_comp" in s_obj else ""
    catch_age_basis = "harvest_rate" in around_pred
    diag_selected = ("selected_numbers" in s_age or "pred[i]" in s_age) and ("n[i] * sel[i]" in s_age or "n[i] * selectivity[i]" in s_age)
    add(level, "2_age_comp_prediction_basis", "inspect",
        f"objective_selected_numbers_basis={pred_selected}; nearby_catch_age_terms={catch_age_basis}; diagnostics_selected_basis={diag_selected}")

    years_pos = s_obj.find("for (std::size_t t = 0; t < years.size(); ++t)")
    obs_pos = s_obj.find("for (const auto &obs", years_pos)
    next_pos = s_obj.find("std::array<T, kAges> next", years_pos)
    rec_pos = s_obj.find("rec_dev", years_pos)
    detail = f"within yearly loop: rec_dev_pos={rec_pos-years_pos if years_pos>=0 and rec_pos>=0 else 'NA'}, obs_loop_pos={obs_pos-years_pos if years_pos>=0 and obs_pos>=0 else 'NA'}, next_update_pos={next_pos-years_pos if years_pos>=0 and next_pos>=0 else 'NA'}"
    if years_pos >= 0 and obs_pos >= 0 and next_pos >= 0:
        detail += "; observations appear before annual survival/update" if obs_pos < next_pos else "; observations appear after annual survival/update"
    add(level, "3_population_update_timing", "inspect", detail)

    has_baranov = all(x in s_obj for x in ["f_a / z_a", "T(1.0) - exp_t(-z_a)", "total_catch_hat"])
    add(level, "4_catch_equation", "pass" if has_baranov else "warn",
        "Baranov-style catch equation found" if has_baranov else "Could not confirm Baranov catch equation")

    has_index = all(x in s_obj for x in ["vulnerable_biomass", "fleet_q * vulnerable_biomass"])
    add(level, "5_index_equation", "pass" if has_index else "warn",
        "index_hat = q * vulnerable biomass pattern found" if has_index else "Could not confirm index equation")

    has_mult = "age_comp_nll" in s_obj and "effective_n" in s_obj and "log_t" in s_obj
    add(level, "6_age_comp_likelihood", "pass" if has_mult else "warn",
        "multinomial-like -N_eff * obs * log(pred) pattern found" if has_mult else "Could not confirm age comp likelihood")

    plus_patterns = [
        "next[last] = next[last] + n[last] * exp_t(-z_last)",
        "n[static_cast<std::size_t>(kAges - 1)] / (T(1.0) - exp_t(-m))",
        "T(1.0) - exp_t(-m_at_age[static_cast<std::size_t>(kAges - 1)])",
    ]
    found_plus = [p for p in plus_patterns if p in s_obj]
    add(level, "7_plus_group", "inspect" if found_plus else "warn",
        f"plus_group_patterns_found={len(found_plus)}; " + " | ".join(found_plus[:2]))

    m_at_age = "m_at_age" in s_obj
    m_scalar = re.search(r'const T m\s*=', s_obj) is not None
    init_uses_m = "exp_t(-m)" in s_obj[:s_obj.find("for (std::size_t t")] if "for (std::size_t t" in s_obj else "exp_t(-m)" in s_obj
    annual_uses_m = "z_a = m +" in s_obj or "z_a = m_at_age" in s_obj
    diag_m_age = "m_at_age" in (s_age + s_init + s_ll + s_ps) or "m_young" in (s_age + s_init + s_ll + s_ps)
    add(level, "8_natural_mortality_consistency", "inspect",
        f"objective_m_at_age={m_at_age}; objective_scalar_m={m_scalar}; init_uses_m={init_uses_m}; annual_uses_m={annual_uses_m}; diagnostics_age_m_aware={diag_m_age}")

with csv_path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["level", "check", "status", "detail"])
    w.writeheader()
    w.writerows(checks)

lines = []
lines.append("Bigeye Objective Operator Audit")
lines.append("===============================")
lines.append("")
lines.append("Purpose")
lines.append("-------")
lines.append("Static audit of BigeyeQuadraObjective::operator() and diagnostics for:")
lines.append("1 parameter layout, 2 age-comp prediction path, 3 population timing,")
lines.append("4 catch equation, 5 index equation, 6 age-comp likelihood,")
lines.append("7 plus group, and 8 natural mortality consistency.")
lines.append("")
for level in levels:
    lines.append(level)
    lines.append("-" * len(level))
    for r in [x for x in checks if x["level"] == level]:
        lines.append(f"{r['check']}: {r['status']}")
        lines.append(f"  {r['detail']}")
    lines.append("")

lines.append("High-risk items to inspect manually")
lines.append("-----------------------------------")
lines.append("A. Age composition is currently predicted from selected numbers at age.")
lines.append("   Verify whether the synthetic observation CSV represents selected survey/fishery")
lines.append("   age composition or catch-at-age composition. If it represents catch-at-age,")
lines.append("   the objective should use N_a * F_a/Z_a * (1-exp(-Z_a)) instead.")
lines.append("")
lines.append("B. Observation timing appears to be before the annual survival/update step.")
lines.append("   Verify whether catches, indices, and age comps are intended to observe beginning-")
lines.append("   of-year population, mid-year population, or catch over the year.")
lines.append("")
lines.append("C. Plus-group initialization and annual update should be checked with M and F.")
lines.append("   Equilibrium plus-group often uses survival into plus group divided by")
lines.append("   1 - survival in plus group. If annual F contributes to survival, make sure")
lines.append("   this is intentional and consistent.")
lines.append("")
lines.append("D. Level 21 age-based M must be consistently applied in equilibrium initialization,")
lines.append("   catch Z, annual survival, plus group, residual diagnostics, and report writers.")
lines.append("")

report.write_text("\n".join(lines) + "\n")
print(f"wrote: {report}")
print(f"wrote: {csv_path}")
PY

cat "$REPORT"
