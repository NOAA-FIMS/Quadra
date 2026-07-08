#!/usr/bin/env bash
set -euo pipefail

# Diagnostic only: no model-code patches.
#
# Purpose:
#   Inspect the current Level 21 fit after the modeling fixes and determine
#   whether remaining non-convergence is dominated by parameter scaling /
#   optimizer line-search behavior rather than a structural model bug.
#
# Usage from repo root:
#   chmod +x diagnose_bigeye_level21_optimizer_scaling.sh
#   ./diagnose_bigeye_level21_optimizer_scaling.sh

ROOT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
OUT="${ROOT}/outputs"
WF="examples/NMFS/pifsc_bigeye_tuna/workflow"
mkdir -p "${WF}"

SUMMARY="${OUT}/bigeye_level21_fit_summary.csv"
GRAD_TXT="${OUT}/bigeye_level21_gradient_by_parameter.txt"
GRAD_CSV="${OUT}/bigeye_level21_gradient_by_parameter.csv"
SANITY="${OUT}/bigeye_level21_parameter_sanity_diagnostics.txt"
RESID="${OUT}/bigeye_level21_age_comp_residual_diagnostics.txt"

REPORT="${WF}/bigeye_level21_optimizer_scaling_diagnostic.txt"
CSV="${WF}/bigeye_level21_optimizer_scaling_diagnostic.csv"

echo "== Ensure current Level 21 outputs exist =="
if [[ ! -f "${SUMMARY}" || ! -f "${GRAD_CSV}" ]]; then
  echo "Level 21 outputs missing; running model first."
  ./run_bigeye_level21_age_based_m_check.sh
fi

echo "== Writing optimizer scaling diagnostic =="
python3 - <<'PY'
from pathlib import Path
import csv
import math
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")
out = root / "outputs"
wf = Path("examples/NMFS/pifsc_bigeye_tuna/workflow")
summary_path = out / "bigeye_level21_fit_summary.csv"
grad_csv_path = out / "bigeye_level21_gradient_by_parameter.csv"
grad_txt_path = out / "bigeye_level21_gradient_by_parameter.txt"
sanity_path = out / "bigeye_level21_parameter_sanity_diagnostics.txt"
resid_path = out / "bigeye_level21_age_comp_residual_diagnostics.txt"

report_path = wf / "bigeye_level21_optimizer_scaling_diagnostic.txt"
csv_path = wf / "bigeye_level21_optimizer_scaling_diagnostic.csv"

def read_key_csv(path):
    d = {}
    if not path.exists():
        return d
    with path.open(newline="") as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                d[row[0]] = row[1]
    return d

summary = read_key_csv(summary_path)

grad_rows = []
if grad_csv_path.exists():
    with grad_csv_path.open(newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            # Expected rows:
            # rank,name,value,fixed_gradient,abs_fixed_gradient
            # or no header with same layout.
            if row[0] == "rank":
                continue
            try:
                rank = int(row[0])
                name = row[1]
                value = float(row[2])
                grad = float(row[3])
                absgrad = float(row[4])
                grad_rows.append((rank, name, value, grad, absgrad))
            except Exception:
                continue

def group_for(name):
    if name in {"log_r0", "log_fbar", "log_q_purse_seine"}:
        return "base_scale"
    if name.startswith("log_m_"):
        return "age_m"
    if name.startswith("init_log_number"):
        return "initial_numbers"
    if name.startswith("logit_sel_longline"):
        return "longline_selectivity"
    if name.startswith("logit_sel_purse_seine"):
        return "purse_seine_selectivity"
    if name.startswith("rec_dev") or "recruit" in name:
        return "recruitment"
    return "other"

groups = {}
for row in grad_rows:
    _, name, value, grad, absgrad = row
    g = group_for(name)
    groups.setdefault(g, []).append(row)

group_summary = []
for g, rows in groups.items():
    abs_sum = sum(r[4] for r in rows)
    max_row = max(rows, key=lambda r: r[4])
    rms = math.sqrt(sum(r[3] ** 2 for r in rows))
    group_summary.append((g, len(rows), abs_sum, rms, max_row[1], max_row[4]))
group_summary.sort(key=lambda x: x[3], reverse=True)

total_rms = math.sqrt(sum(r[3] ** 2 for r in grad_rows)) if grad_rows else float("nan")
base_rms = next((x[3] for x in group_summary if x[0] == "base_scale"), 0.0)
base_share = base_rms / total_rms if total_rms and math.isfinite(total_rms) else float("nan")

top = sorted(grad_rows, key=lambda r: r[4], reverse=True)[:25]

# Parse stress/residual snippets.
sanity_text = sanity_path.read_text(errors="replace") if sanity_path.exists() else ""
resid_text = resid_path.read_text(errors="replace") if resid_path.exists() else ""

prior_block = ""
m = re.search(r"Prior penalty by block\n-+\n(.*?)(?:\n\n|\Z)", sanity_text, re.S)
if m:
    prior_block = m.group(1).strip()

stress_block = ""
m = re.search(r"Initial-number stress summary\n-+\n(.*?)(?:\n\n|\Z)", sanity_text, re.S)
if m:
    stress_block = m.group(1).strip()

resid_block = ""
m = re.search(r"Fleet summary\n-+\n(.*?)(?:\n\n|\Z)", resid_text, re.S)
if m:
    resid_block = m.group(1).strip()

objective = summary.get("objective", "NA")
grad_norm = summary.get("grad_norm", "NA")
converged = summary.get("converged", "NA")

interpretation = []
try:
    gn = float(grad_norm)
except Exception:
    gn = float("nan")

if math.isfinite(gn) and gn < 0.05:
    interpretation.append("The final gradient is small but above the requested tolerance; this is consistent with optimizer/line-search polishing trouble, not necessarily a broken objective.")
if base_share > 0.50:
    interpretation.append("More than half of the gradient RMS is in base-scale parameters; prioritize scale-aware polishing of log_r0, log_fbar, and log_q before changing biology.")
elif base_share > 0.25:
    interpretation.append("Base-scale parameters are a meaningful part of the remaining gradient; a scale-aware polish is still worth testing.")
else:
    interpretation.append("The remaining gradient is not dominated by base-scale parameters; inspect the leading age/selectivity block before optimizer work.")

if "age_based_m_prior_nll" in prior_block and "initial_numbers_prior_nll" in prior_block:
    interpretation.append("The prior block gives context for whether the fit is paying penalty in M, initial numbers, or selectivity. Compare these against objective changes before patching model equations.")

recommendations = [
    "Run a restart from the final iterate with a looser line-search epsilon or a small gradient-descent/Barzilai-Borwein polish.",
    "Try parameter scaling: base log parameters can be stepped separately from logits and initial-number deviations.",
    "Do not patch plus-group or M again based only on this run; fixed-M was worse and plus dynamics looked stable.",
    "Next high-value check: add a one-step directional descent diagnostic, reporting f(theta - alpha*g) for alpha in {1e-4, 3e-4, 1e-3, 3e-3, 1e-2}.",
]

with report_path.open("w") as f:
    f.write("Bigeye Level 21 Optimizer Scaling Diagnostic\n")
    f.write("===========================================\n\n")
    f.write("Fit summary\n")
    f.write("-----------\n")
    f.write(f"objective,{objective}\n")
    f.write(f"grad_norm,{grad_norm}\n")
    f.write(f"converged,{converged}\n\n")

    f.write("Gradient group summary\n")
    f.write("----------------------\n")
    f.write("group,n_params,abs_gradient_sum,rms_gradient,max_gradient_param,max_abs_gradient\n")
    for row in group_summary:
        f.write(",".join([row[0], str(row[1]), f"{row[2]:.12g}", f"{row[3]:.12g}", row[4], f"{row[5]:.12g}"]) + "\n")
    f.write(f"\ntotal_gradient_rms,{total_rms:.12g}\n")
    f.write(f"base_scale_rms_share,{base_share:.12g}\n\n")

    f.write("Top gradients\n")
    f.write("-------------\n")
    f.write("rank,name,value,fixed_gradient,abs_fixed_gradient,group\n")
    for rank, name, value, grad, absgrad in top:
        f.write(f"{rank},{name},{value:.12g},{grad:.12g},{absgrad:.12g},{group_for(name)}\n")

    if prior_block:
        f.write("\nPrior penalty by block\n")
        f.write("----------------------\n")
        f.write(prior_block + "\n")
    if stress_block:
        f.write("\nInitial-number stress summary\n")
        f.write("-----------------------------\n")
        f.write(stress_block + "\n")
    if resid_block:
        f.write("\nResidual fleet summary\n")
        f.write("----------------------\n")
        f.write(resid_block + "\n")

    f.write("\nInterpretation\n")
    f.write("--------------\n")
    for s in interpretation:
        f.write(f"- {s}\n")

    f.write("\nRecommended next diagnostic\n")
    f.write("---------------------------\n")
    for s in recommendations:
        f.write(f"- {s}\n")

with csv_path.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["section", "group", "n_params", "abs_gradient_sum", "rms_gradient", "max_gradient_param", "max_abs_gradient"])
    for row in group_summary:
        w.writerow(["gradient_group", row[0], row[1], row[2], row[3], row[4], row[5]])
    w.writerow(["summary", "total_gradient_rms", "", "", total_rms, "", ""])
    w.writerow(["summary", "base_scale_rms_share", "", "", base_share, "", ""])
    w.writerow([])
    w.writerow(["section", "rank", "name", "value", "fixed_gradient", "abs_fixed_gradient", "group"])
    for rank, name, value, grad, absgrad in top:
        w.writerow(["top_gradient", rank, name, value, grad, absgrad, group_for(name)])

print(f"wrote: {report_path}")
print(f"wrote: {csv_path}")
PY

echo
echo "== Preview =="
sed -n '1,120p' "${REPORT}"

echo
echo "Done."
echo "Report:"
echo "  ${REPORT}"
echo "CSV:"
echo "  ${CSV}"
