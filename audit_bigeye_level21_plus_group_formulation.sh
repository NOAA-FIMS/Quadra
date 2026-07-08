#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna"
L21="$ROOT/level21_age_based_natural_mortality_diagnostic"
OUT="$L21/outputs"

if [[ ! -f "$OUT/bigeye_level21_plus_group_audit.csv" ]]; then
  echo "plus-group audit CSV not found; running Level 21 check first"
  ./run_bigeye_level21_age_based_m_check.sh
fi

python3 - <<'PY'
from pathlib import Path
import csv, math, statistics

ROOT = Path("examples/NMFS/pifsc_bigeye_tuna")
L21 = ROOT / "level21_age_based_natural_mortality_diagnostic"
out = L21 / "outputs"
csv_in = out / "bigeye_level21_plus_group_audit.csv"
txt_out = out / "bigeye_level21_plus_group_formulation_audit.txt"
csv_out = out / "bigeye_level21_plus_group_formulation_audit.csv"

initial = []
plus = []
with csv_in.open() as f:
    r = csv.DictReader(f)
    for row in r:
        if row["section"] == "initial":
            initial.append(row)
        elif row["section"] == "plus_dynamics":
            plus.append(row)

def flt(row, key):
    try:
        return float(row[key])
    except Exception:
        return float("nan")

init_by_age = {int(float(r["age"])): r for r in initial if r["age"] != "NA"}
age9 = init_by_age.get(9)
age10 = init_by_age.get(10)
if not age9 or not age10:
    raise SystemExit("Could not find initial age 9/10 rows in plus-group audit CSV")

n9 = flt(age9, "n")
n10_current = flt(age10, "n")
m10 = flt(age10, "m")
z10 = flt(age10, "z")
ll10 = flt(age10, "ll_sel")
ps10 = flt(age10, "ps_sel")

s_m10 = math.exp(-m10)
s_z10 = math.exp(-z10)
n10_equil_m = n9 * s_m10 / max(1.0 - s_m10, 1e-12)
n10_equil_z = n9 * s_z10 / max(1.0 - s_z10, 1e-12)

plus_survival = [flt(r, "plus_survival") for r in plus]
plus_inflow = [flt(r, "plus_inflow") for r in plus]
plus_next = [flt(r, "plus_next") for r in plus]
plus_n = [flt(r, "n") for r in plus]
ratio_next_current = [pn/n if n and math.isfinite(pn) and math.isfinite(n) else float("nan") for pn,n in zip(plus_next, plus_n)]
ratio_inflow_surv = [i/s if s and math.isfinite(i) and math.isfinite(s) else float("nan") for i,s in zip(plus_inflow, plus_survival)]

def clean(xs):
    return [x for x in xs if math.isfinite(x)]
def mean(xs):
    xs=clean(xs)
    return statistics.mean(xs) if xs else float("nan")
def mn(xs):
    xs=clean(xs)
    return min(xs) if xs else float("nan")
def mx(xs):
    xs=clean(xs)
    return max(xs) if xs else float("nan")

rows = [
    ["metric","value","interpretation"],
    ["n9_initial", n9, "age-9 abundance entering terminal age"],
    ["n10_current_initial", n10_current, "current terminal abundance after plus initialization and init dev"],
    ["m10", m10, "terminal natural mortality"],
    ["z10", z10, "terminal total mortality using total fleet selectivity"],
    ["ll_sel10", ll10, "terminal longline selectivity"],
    ["ps_sel10", ps10, "terminal purse-seine selectivity"],
    ["survival_m10", s_m10, "M-only terminal survival"],
    ["survival_z10", s_z10, "M+F terminal survival"],
    ["equil_n10_m_only", n10_equil_m, "N9*exp(-M10)/(1-exp(-M10))"],
    ["equil_n10_z_total", n10_equil_z, "N9*exp(-Z10)/(1-exp(-Z10))"],
    ["current_minus_m_equil", n10_current - n10_equil_m, "current terminal abundance minus M-only equilibrium"],
    ["current_minus_z_equil", n10_current - n10_equil_z, "current terminal abundance minus Z-total equilibrium"],
    ["z_equil_minus_m_equil", n10_equil_z - n10_equil_m, "effect of including fishing mortality in plus equilibrium"],
    ["mean_annual_plus_inflow", mean(plus_inflow), "age-9 inflow into plus group"],
    ["mean_annual_plus_survival", mean(plus_survival), "surviving plus-group fish"],
    ["mean_annual_plus_next", mean(plus_next), "next-year terminal abundance"],
    ["mean_plus_next_over_current", mean(ratio_next_current), "near 1 means terminal age is dynamically stable"],
    ["min_plus_next_over_current", mn(ratio_next_current), ""],
    ["max_plus_next_over_current", mx(ratio_next_current), ""],
    ["mean_inflow_over_plus_survival", mean(ratio_inflow_surv), "large values mean inflow dominates; small values mean resident plus dominates"],
]

with csv_out.open("w", newline="") as f:
    csv.writer(f).writerows(rows)

with txt_out.open("w") as f:
    f.write("Level 21 Plus-Group Formulation Audit\n")
    f.write("=====================================\n\n")
    f.write("Key formulas checked\n--------------------\n")
    f.write("M-only equilibrium:     Nplus = N9 * exp(-M10) / (1 - exp(-M10))\n")
    f.write("Total-Z equilibrium:    Nplus = N9 * exp(-Z10) / (1 - exp(-Z10))\n")
    f.write("Annual update observed: Nplus_next = N9 * exp(-Z9) + Nplus * exp(-Z10)\n\n")
    f.write("Results\n-------\n")
    for k,v,interp in rows[1:]:
        f.write(f"{k},{v},{interp}\n")
    f.write("\nInterpretation\n--------------\n")
    f.write("If current_initial is close to M-only equilibrium but far from Z-total equilibrium,\n")
    f.write("then initialization and annual transition are not using the same mortality basis.\n")
    f.write("That is not always wrong for an unfished initial equilibrium, but it can create a\n")
    f.write("transient in the terminal age once fishing starts. If the model starts in an exploited\n")
    f.write("period or the synthetic data were generated with fishing active at t0, Z-total may be\n")
    f.write("the more internally consistent initialization.\n")

print(f"wrote: {txt_out}")
print(f"wrote: {csv_out}")
print(txt_out.read_text())
PY
