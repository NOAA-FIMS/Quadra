#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna"
L20="$ROOT/level20_longline_selectivity_regularization_scan"
L21="$ROOT/level21_age_based_natural_mortality_diagnostic"
OUT="$ROOT/workflow/bigeye_level20_vs_level21_gradient_comparison"

if [[ ! -d "$L20" ]]; then
  echo "ERROR: missing $L20"
  exit 1
fi
if [[ ! -d "$L21" ]]; then
  echo "ERROR: missing $L21"
  exit 1
fi

echo "== Ensure Level 21 gradient-by-parameter output exists =="
if [[ ! -f "$L21/outputs/bigeye_level21_gradient_by_parameter.csv" ]]; then
  if [[ -x ./run_bigeye_level21_gradient_by_parameter_check.sh ]]; then
    ./run_bigeye_level21_gradient_by_parameter_check.sh
  else
    echo "ERROR: missing Level 21 gradient CSV and run_bigeye_level21_gradient_by_parameter_check.sh is not executable"
    exit 1
  fi
fi

echo "== Try to generate Level 20 gradient-by-parameter output if absent =="
if [[ ! -f "$L20/outputs/bigeye_level20_gradient_by_parameter.csv" ]]; then
  echo "Level 20 gradient-by-parameter CSV not found."
  echo "Creating a lightweight comparison from Level 20 parameter sanity + Level 21 gradients only."
fi

python3 - <<'PY'
from pathlib import Path
import csv
import math

ROOT = Path("examples/NMFS/pifsc_bigeye_tuna")
L20 = ROOT / "level20_longline_selectivity_regularization_scan"
L21 = ROOT / "level21_age_based_natural_mortality_diagnostic"
OUT = ROOT / "workflow"
OUT.mkdir(parents=True, exist_ok=True)

txt_path = OUT / "bigeye_level20_vs_level21_gradient_comparison.txt"
csv_path = OUT / "bigeye_level20_vs_level21_gradient_comparison.csv"

def read_key_value_csv(path):
    d = {}
    if not path.exists():
        return d
    with path.open() as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                d[row[0]] = row[1]
    return d

def read_gradient_csv(path):
    rows = []
    if not path.exists():
        return rows
    with path.open() as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                row["value_num"] = float(row.get("value", "nan"))
            except Exception:
                row["value_num"] = math.nan
            try:
                row["grad_num"] = float(row.get("fixed_gradient", "nan"))
            except Exception:
                row["grad_num"] = math.nan
            try:
                row["abs_grad_num"] = float(row.get("abs_fixed_gradient", "nan"))
            except Exception:
                row["abs_grad_num"] = abs(row["grad_num"]) if math.isfinite(row["grad_num"]) else math.nan
            rows.append(row)
    return rows

l20_fit = read_key_value_csv(L20 / "outputs/bigeye_level20_fit_summary.csv")
l21_fit = read_key_value_csv(L21 / "outputs/bigeye_level21_fit_summary.csv")
l21_grad = read_gradient_csv(L21 / "outputs/bigeye_level21_gradient_by_parameter.csv")
l20_grad = read_gradient_csv(L20 / "outputs/bigeye_level20_gradient_by_parameter.csv")

def f(d, k):
    try:
        return float(d.get(k, "nan"))
    except Exception:
        return math.nan

summary_rows = [
    ["metric", "level20", "level21", "delta_21_minus_20"],
    ["objective", l20_fit.get("objective",""), l21_fit.get("objective",""), f(l21_fit,"objective") - f(l20_fit,"objective")],
    ["grad_norm", l20_fit.get("grad_norm",""), l21_fit.get("grad_norm",""), f(l21_fit,"grad_norm") - f(l20_fit,"grad_norm")],
    ["converged", l20_fit.get("converged",""), l21_fit.get("converged",""), ""],
    ["log_r0", l20_fit.get("log_r0",""), l21_fit.get("log_r0",""), f(l21_fit,"log_r0") - f(l20_fit,"log_r0")],
    ["r0", l20_fit.get("r0",""), l21_fit.get("r0",""), f(l21_fit,"r0") - f(l20_fit,"r0")],
    ["log_fbar", l20_fit.get("log_fbar",""), l21_fit.get("log_fbar",""), f(l21_fit,"log_fbar") - f(l20_fit,"log_fbar")],
    ["fbar", l20_fit.get("fbar",""), l21_fit.get("fbar",""), f(l21_fit,"fbar") - f(l20_fit,"fbar")],
    ["log_q_purse_seine", l20_fit.get("log_q_purse_seine",""), l21_fit.get("log_q_purse_seine",""), f(l21_fit,"log_q_purse_seine") - f(l20_fit,"log_q_purse_seine")],
    ["q_purse_seine", l20_fit.get("q_purse_seine",""), l21_fit.get("q_purse_seine",""), f(l21_fit,"q_purse_seine") - f(l20_fit,"q_purse_seine")],
    ["log_m_young_offset", "", l21_fit.get("log_m_young_offset",""), ""],
    ["log_m_old_offset", "", l21_fit.get("log_m_old_offset",""), ""],
]

# Add init multipliers
for a in range(1, 11):
    k = f"init_number_multiplier_age_{a}"
    summary_rows.append([k, l20_fit.get(k,""), l21_fit.get(k,""), f(l21_fit,k)-f(l20_fit,k)])

# If Level 20 gradient exists, compare by name. Otherwise report Level 21 top gradients only.
comparison_rows = []
if l20_grad:
    l20_by_name = {r.get("name",""): r for r in l20_grad}
    l21_by_name = {r.get("name",""): r for r in l21_grad}
    names = sorted(set(l20_by_name) | set(l21_by_name))
    for name in names:
        g20 = l20_by_name.get(name, {})
        g21 = l21_by_name.get(name, {})
        v20 = g20.get("value_num", math.nan)
        v21 = g21.get("value_num", math.nan)
        a20 = g20.get("abs_grad_num", math.nan)
        a21 = g21.get("abs_grad_num", math.nan)
        comparison_rows.append({
            "name": name,
            "level20_value": v20,
            "level21_value": v21,
            "delta_value": v21 - v20 if math.isfinite(v20) and math.isfinite(v21) else "",
            "level20_abs_gradient": a20,
            "level21_abs_gradient": a21,
            "delta_abs_gradient": a21 - a20 if math.isfinite(a20) and math.isfinite(a21) else "",
        })
    comparison_rows.sort(key=lambda r: (float(r["level21_abs_gradient"]) if isinstance(r["level21_abs_gradient"], float) and math.isfinite(r["level21_abs_gradient"]) else -1), reverse=True)
else:
    for r in sorted(l21_grad, key=lambda x: x.get("abs_grad_num", math.nan), reverse=True):
        comparison_rows.append({
            "name": r.get("name",""),
            "level20_value": "",
            "level21_value": r.get("value",""),
            "delta_value": "",
            "level20_abs_gradient": "",
            "level21_abs_gradient": r.get("abs_fixed_gradient",""),
            "delta_abs_gradient": "",
        })

with csv_path.open("w", newline="") as fcsv:
    w = csv.writer(fcsv)
    w.writerow(["section","name","level20","level21","delta_21_minus_20","level20_abs_gradient","level21_abs_gradient","delta_abs_gradient"])
    for row in summary_rows[1:]:
        w.writerow(["summary", row[0], row[1], row[2], row[3], "", "", ""])
    for r in comparison_rows:
        w.writerow(["gradient", r["name"], r["level20_value"], r["level21_value"], r["delta_value"], r["level20_abs_gradient"], r["level21_abs_gradient"], r["delta_abs_gradient"]])

top = comparison_rows[:20]

with txt_path.open("w") as txt:
    txt.write("Bigeye Level 20 vs Level 21 Gradient Comparison\n")
    txt.write("================================================\n\n")
    txt.write("Purpose\n-------\n")
    txt.write("Compare Level 20 and Level 21 fixed-effect behavior to determine whether\n")
    txt.write("the high Level 21 gradient is caused by age-based M coupling, initial-number\n")
    txt.write("stress, or a remaining parameter-layout/reporting mismatch.\n\n")

    txt.write("Fit summary comparison\n----------------------\n")
    for row in summary_rows[1:]:
        txt.write(f"{row[0]},{row[1]},{row[2]},{row[3]}\n")

    txt.write("\nTop Level 21 gradient rows\n--------------------------\n")
    txt.write("name,level20_value,level21_value,delta_value,level20_abs_gradient,level21_abs_gradient,delta_abs_gradient\n")
    for r in top:
        txt.write(f"{r['name']},{r['level20_value']},{r['level21_value']},{r['delta_value']},{r['level20_abs_gradient']},{r['level21_abs_gradient']},{r['delta_abs_gradient']}\n")

    txt.write("\nInterpretation\n--------------\n")
    if not l20_grad:
        txt.write("Level 20 gradient-by-parameter CSV was not found, so this report ranks Level 21 gradients only.\n")
        txt.write("The next useful step is to add the same gradient-by-parameter writer to Level 20, then compare matched parameter names.\n")
    else:
        txt.write("Rows with large positive delta_abs_gradient identify parameters that became harder to optimize after adding age-based M.\n")
    txt.write("If top gradients concentrate in initial numbers and M offsets while residuals remain good, the model likely has a near-equivalent ridge between initial numbers, age-specific M, and selectivity.\n")

print(f"wrote: {txt_path}")
print(f"wrote: {csv_path}")

print("\n== Preview ==")
print(txt_path.read_text().split("Top Level 21 gradient rows")[0])
print("Top Level 21 gradient rows")
print("--------------------------")
for line in txt_path.read_text().split("Top Level 21 gradient rows\n--------------------------\n",1)[1].splitlines()[:22]:
    print(line)
PY
