#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
OUT="${ROOT}/outputs"
WF="examples/NMFS/pifsc_bigeye_tuna/workflow"
mkdir -p "${WF}"

SUMMARY="${OUT}/bigeye_level21_fit_summary.csv"
GRAD_CSV="${OUT}/bigeye_level21_gradient_by_parameter.csv"
REPORT="${WF}/bigeye_level21_directional_descent_diagnostic.txt"
CSV="${WF}/bigeye_level21_directional_descent_diagnostic.csv"

echo "== Ensure current Level 21 outputs exist =="
if [[ ! -f "${SUMMARY}" || ! -f "${GRAD_CSV}" ]]; then
  ./run_bigeye_level21_age_based_m_check.sh
fi

echo "== Build first-order descent table =="
python3 - <<'PY'
from pathlib import Path
import csv, math

out = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs")
wf = Path("examples/NMFS/pifsc_bigeye_tuna/workflow")
summary_path = out / "bigeye_level21_fit_summary.csv"
grad_path = out / "bigeye_level21_gradient_by_parameter.csv"
report_path = wf / "bigeye_level21_directional_descent_diagnostic.txt"
csv_path = wf / "bigeye_level21_directional_descent_diagnostic.csv"

def read_key_csv(path):
    d = {}
    with path.open(newline="") as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                d[row[0]] = row[1]
    return d

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

summary = read_key_csv(summary_path)
rows = []
with grad_path.open(newline="") as f:
    for row in csv.reader(f):
        if not row or row[0] == "rank":
            continue
        try:
            rows.append({
                "rank": int(row[0]),
                "name": row[1],
                "value": float(row[2]),
                "grad": float(row[3]),
                "absgrad": float(row[4]),
                "group": group_for(row[1]),
            })
        except Exception:
            pass

alphas = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3]
groups = ["all", "age_m", "base_scale", "purse_seine_selectivity",
          "longline_selectivity", "initial_numbers"]

diagnostics = []
for gname in groups:
    use = rows if gname == "all" else [r for r in rows if r["group"] == gname]
    if not use:
        continue
    grad_sq = sum(r["grad"] ** 2 for r in use)
    grad_rms = math.sqrt(grad_sq)
    max_row = max(use, key=lambda r: abs(r["grad"]))
    for alpha in alphas:
        diagnostics.append({
            "direction": gname,
            "alpha": alpha,
            "n_params": len(use),
            "grad_rms": grad_rms,
            "grad_sq": grad_sq,
            "predicted_delta": -alpha * grad_sq,
            "max_abs_step": alpha * abs(max_row["grad"]),
            "max_grad_param": max_row["name"],
            "max_abs_grad": abs(max_row["grad"]),
        })

recommended = []
for gname in groups:
    cand = [d for d in diagnostics if d["direction"] == gname and d["max_abs_step"] <= 1e-5]
    if not cand:
        cand = [d for d in diagnostics if d["direction"] == gname]
    if cand:
        recommended.append(max(cand, key=lambda d: abs(d["predicted_delta"])))

with report_path.open("w") as f:
    f.write("Bigeye Level 21 Directional Descent Diagnostic\n")
    f.write("=============================================\n\n")
    f.write("Fit summary\n-----------\n")
    f.write(f"objective,{summary.get('objective','NA')}\n")
    f.write(f"grad_norm,{summary.get('grad_norm','NA')}\n")
    f.write(f"converged,{summary.get('converged','NA')}\n\n")
    f.write("Important limitation\n--------------------\n")
    f.write("This script computes first-order expected decreases from reported gradients.\n")
    f.write("The actual f(theta - alpha*g) check requires a tiny C++ objective-evaluation driver.\n\n")
    f.write("Recommended actual objective checks\n-----------------------------------\n")
    f.write("direction,alpha,max_abs_step,predicted_first_order_delta,max_grad_param,max_abs_grad\n")
    for d in recommended:
        f.write(f"{d['direction']},{d['alpha']:.12g},{d['max_abs_step']:.12g},"
                f"{d['predicted_delta']:.12g},{d['max_grad_param']},{d['max_abs_grad']:.12g}\n")
    f.write("\nFull first-order table\n----------------------\n")
    f.write("direction,alpha,n_params,grad_rms,grad_sq,predicted_first_order_delta,max_abs_step,max_grad_param,max_abs_grad\n")
    for d in diagnostics:
        f.write(f"{d['direction']},{d['alpha']:.12g},{d['n_params']},{d['grad_rms']:.12g},"
                f"{d['grad_sq']:.12g},{d['predicted_delta']:.12g},{d['max_abs_step']:.12g},"
                f"{d['max_grad_param']},{d['max_abs_grad']:.12g}\n")
    f.write("\nInterpretation\n--------------\n")
    f.write("- If actual f(theta - alpha*g) decreases for small alpha, the gradient is usable and L-BFGS needs polishing/scaling.\n")
    f.write("- If actual f(theta - alpha*g) does not decrease even for tiny alpha, suspect stale parameter/gradient output or gradient-objective mismatch.\n")
    f.write("- Test full, age_m-only, and base_scale-only directions first.\n")

with csv_path.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["direction","alpha","n_params","grad_rms","grad_sq","predicted_first_order_delta","max_abs_step","max_grad_param","max_abs_grad"])
    for d in diagnostics:
        w.writerow([d["direction"], d["alpha"], d["n_params"], d["grad_rms"], d["grad_sq"],
                    d["predicted_delta"], d["max_abs_step"], d["max_grad_param"], d["max_abs_grad"]])

print(f"wrote: {report_path}")
print(f"wrote: {csv_path}")
PY

echo
echo "== Preview =="
sed -n '1,120p' "${REPORT}"
