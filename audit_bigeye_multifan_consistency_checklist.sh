#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna"
L20="${BASE}/level20_longline_selectivity_regularization_scan"
L21="${BASE}/level21_age_based_natural_mortality_diagnostic"
L23="${BASE}/level23_longline_selectivity_smoothness_scan"
WF="${BASE}/workflow"

mkdir -p "${WF}"

REPORT="${WF}/bigeye_multifan_consistency_checklist.txt"
CSV="${WF}/bigeye_multifan_consistency_checklist.csv"

python3 - <<'PY'
from pathlib import Path
import csv
import re

base = Path("examples/NMFS/pifsc_bigeye_tuna")
levels = {
    "level20": base / "level20_longline_selectivity_regularization_scan",
    "level21": base / "level21_age_based_natural_mortality_diagnostic",
    "level23": base / "level23_longline_selectivity_smoothness_scan",
}
wf = base / "workflow"
report = wf / "bigeye_multifan_consistency_checklist.txt"
csvout = wf / "bigeye_multifan_consistency_checklist.csv"

def read_text(path):
    return path.read_text(errors="replace") if path.exists() else ""

def read_key_csv(path):
    d = {}
    if path.exists():
        with path.open(newline="") as f:
            for row in csv.reader(f):
                if len(row) >= 2:
                    d[row[0]] = row[1]
    return d

def block_extract(text, heading):
    m = re.search(re.escape(heading) + r"\n-+\n(.*?)(?:\n\n|\Z)", text, re.S)
    return m.group(1).strip() if m else ""

def summarize_level(name, root):
    obj = read_text(root / "objective" / "bigeye_quadra_objective.hpp")
    summary = read_key_csv(root / "outputs" / f"bigeye_{name}_fit_summary.csv")

    has_age_m_params = "log_m_young_offset" in obj or "log_m_old_offset" in obj
    has_ll_template_prior = "ll_template" in obj and ("sigma_ll_sel_dev" in obj or "BIGEYE_LL_SEL_SIGMA" in obj)
    has_ps_template_prior = "ps_template" in obj and "sigma_ps_sel_dev" in obj
    has_total_f = "f_total_a" in obj or ("f_longline_a" in obj and "f_purse_seine_a" in obj)
    has_fleet_baranov = "catch_at_age_longline" in obj and "catch_at_age_purse_seine" in obj
    has_old_catch_share = "fleet_catch_share" in obj or "total_catch_hat" in obj
    comp_uses_catch_at_age = "fleet_catch_at_age" in obj and "pred_age_comp[i] = fleet_catch_at_age[i]" in obj
    comp_uses_vulnerable_numbers = re.search(r"pred_age_comp\[i\]\s*=\s*n\[i\]\s*\*\s*sel\[i\]", obj) is not None
    plus_group = "next[last] = next[last] + n[last]" in obj
    plus_init_m_only = "1.0) - exp_t(-m_at_age" in obj or "T(1.0) - exp_t(-m_at_age" in obj

    return {
        "level": name,
        "objective": summary.get("objective", ""),
        "grad_norm": summary.get("grad_norm", ""),
        "converged": summary.get("converged", ""),
        "m_convention": "estimated_age_m" if has_age_m_params else "fixed_m_like",
        "selectivity_regularization": "template_prior" if (has_ll_template_prior or has_ps_template_prior) else "unclear",
        "fleet_specific_baranov": "yes" if has_fleet_baranov else "no",
        "total_f_in_survival": "yes" if has_total_f else "unclear/no",
        "old_catch_share_leftover": "yes" if has_old_catch_share else "no",
        "composition_basis": "catch_at_age" if comp_uses_catch_at_age else ("vulnerable_numbers" if comp_uses_vulnerable_numbers else "unclear"),
        "plus_group_dynamics": "plus_accumulator" if plus_group else "unclear",
        "plus_initialization": "m_only_equilibrium" if plus_init_m_only else "unclear",
    }

rows = [summarize_level(k, v) for k, v in levels.items()]

l21 = levels["level21"]
l21_sanity = read_text(l21 / "outputs" / "bigeye_level21_parameter_sanity_diagnostics.txt")
l21_grad = read_text(l21 / "outputs" / "bigeye_level21_gradient_by_parameter.txt")
l21_scale = read_text(wf / "bigeye_level21_optimizer_scaling_diagnostic.txt")

prior_block = block_extract(l21_sanity, "Prior penalty by block")
init_block = block_extract(l21_sanity, "Initial-number stress summary")
grad_group_block = block_extract(l21_scale, "Gradient group summary")
top_grad_block = block_extract(l21_grad, "Top gradients")

assessment = [
    ("M-at-age", "Published-style convention used fixed M-at-age.", "Level 21 estimates young/old M offsets; fixed-M diagnostic was worse, so M helps fit but may be confounded."),
    ("Selectivity", "Published-style convention used smoothed/constrained selectivity.", "Template priors help, but Level 21 still pays large selectivity penalties; Level 23 smoothness scan is the right comparator."),
    ("Catch equation", "Fleet-specific Baranov catch is preferred.", "Level 21 appears to use fleet-specific catch-at-age and no old fleet_catch_share/total_catch_hat path."),
    ("Composition likelihood", "Fishery age comps should usually use catch-at-age basis.", "Level 21 appears to use fleet catch-at-age for predicted compositions."),
    ("Plus group", "Plus group receives age-9 inflow plus surviving terminal fish.", "Level 21 appears to use accumulator dynamics; M-only initialization is defensible for unfished initial state."),
    ("Optimizer", "Check actual f(theta - alpha*g).", "Remaining gradient has age-M/base-scale/selectivity components; directional descent is the decisive next test."),
]

with report.open("w") as f:
    f.write("Bigeye MULTIFAN/A-SCALA Consistency Checklist\n")
    f.write("============================================\n\n")
    f.write("Level summary\n-------------\n")
    f.write("level,objective,grad_norm,converged,m_convention,selectivity_regularization,fleet_specific_baranov,total_f_in_survival,old_catch_share_leftover,composition_basis,plus_group_dynamics,plus_initialization\n")
    for r in rows:
        f.write(",".join(str(r[k]) for k in [
            "level","objective","grad_norm","converged","m_convention",
            "selectivity_regularization","fleet_specific_baranov","total_f_in_survival",
            "old_catch_share_leftover","composition_basis","plus_group_dynamics",
            "plus_initialization"
        ]) + "\n")

    f.write("\nInterpretive checklist\n----------------------\n")
    f.write("topic,published_or_target_convention,current_interpretation\n")
    for topic, convention, interp in assessment:
        f.write(f"{topic},{convention},{interp}\n")

    if prior_block:
        f.write("\nLevel 21 prior penalty by block\n-------------------------------\n")
        f.write(prior_block + "\n")
    if init_block:
        f.write("\nLevel 21 initial-number stress\n------------------------------\n")
        f.write(init_block + "\n")
    if grad_group_block:
        f.write("\nLevel 21 gradient group summary\n-------------------------------\n")
        f.write(grad_group_block + "\n")
    if top_grad_block:
        f.write("\nLevel 21 top gradients\n----------------------\n")
        f.write("\n".join(top_grad_block.splitlines()[:15]) + "\n")

    f.write("\nRecommended next actions\n------------------------\n")
    f.write("1. Run actual directional descent: f(theta - alpha*g) for full, age_m-only, and base-scale-only directions.\n")
    f.write("2. Compare Level 21 against Level 23 smoothness scan before changing M or plus-group again.\n")
    f.write("3. Treat age-M as useful but potentially non-identifiable unless directional descent/profile diagnostics support it.\n")
    f.write("4. Keep fleet-specific Baranov and catch-at-age composition as the preferred current formulation.\n")

with csvout.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["section","level","field","value"])
    for r in rows:
        for k, v in r.items():
            if k != "level":
                w.writerow(["level_summary", r["level"], k, v])
    w.writerow([])
    w.writerow(["section","topic","published_or_target_convention","current_interpretation"])
    for topic, convention, interp in assessment:
        w.writerow(["checklist", topic, convention, interp])

print(f"wrote: {report}")
print(f"wrote: {csvout}")
PY

echo "== Preview =="
sed -n '1,180p' "${REPORT}"
