#!/usr/bin/env bash
set -euo pipefail

OUT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs"
A="$OUT/bigeye_level23_longline_prediction_decomposition_smooth_0.csv"
B="$OUT/bigeye_level23_longline_prediction_decomposition_smooth_0.01.csv"
REPORT="$OUT/bigeye_level23_smooth_0_vs_001_comparison.txt"
CSV="$OUT/bigeye_level23_smooth_0_vs_001_comparison.csv"

if [[ ! -f "$A" || ! -f "$B" ]]; then
  echo "ERROR: missing decomposition files. Re-run:"
  echo "  ./run_bigeye_level23_longline_selectivity_smoothness_scan.sh"
  exit 1
fi

python3 - <<'PY'
from pathlib import Path
import csv
from collections import defaultdict

out = Path("examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs")
a_path = out / "bigeye_level23_longline_prediction_decomposition_smooth_0.csv"
b_path = out / "bigeye_level23_longline_prediction_decomposition_smooth_0.01.csv"
txt_path = out / "bigeye_level23_smooth_0_vs_001_comparison.txt"
csv_path = out / "bigeye_level23_smooth_0_vs_001_comparison.csv"

def read_rows(path):
    rows = []
    with path.open() as f:
        for row in csv.DictReader(f):
            if not row or row.get("year") in (None, "", "year"):
                continue
            try:
                sel = row.get("longline_selectivity", row.get("selectivity", "nan"))
                rows.append({
                    "year": int(float(row["year"])),
                    "age": int(float(row["age"])),
                    "observed": float(row["observed_comp"]),
                    "predicted": float(row["predicted_comp"]),
                    "residual": float(row["residual"]),
                    "abs_residual": float(row["abs_residual"]),
                    "selectivity": float(sel),
                    "n_at_age": float(row["n_at_age"]),
                    "selected_numbers": float(row["selected_numbers"]),
                })
            except Exception:
                pass
    return rows

a = read_rows(a_path)
b = read_rows(b_path)

def summarize(rows):
    by_age = defaultdict(list)
    for r in rows:
        by_age[r["age"]].append(r)
    out_rows = []
    for age in sorted(by_age):
        vals = by_age[age]
        out_rows.append({
            "age": age,
            "mean_obs": sum(v["observed"] for v in vals)/len(vals),
            "mean_pred": sum(v["predicted"] for v in vals)/len(vals),
            "mean_abs": sum(v["abs_residual"] for v in vals)/len(vals),
            "max_abs": max(v["abs_residual"] for v in vals),
            "mean_sel": sum(v["selectivity"] for v in vals)/len(vals),
        })
    return out_rows

sa = summarize(a)
sb = summarize(b)
b_by_age = {r["age"]: r for r in sb}

with csv_path.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "age",
        "mean_abs_lambda_0",
        "mean_abs_lambda_0.01",
        "delta_mean_abs_0.01_minus_0",
        "max_abs_lambda_0",
        "max_abs_lambda_0.01",
        "mean_selectivity_lambda_0",
        "mean_selectivity_lambda_0.01",
        "delta_selectivity_0.01_minus_0",
    ])
    for r in sa:
        q = b_by_age[r["age"]]
        w.writerow([
            r["age"],
            r["mean_abs"],
            q["mean_abs"],
            q["mean_abs"] - r["mean_abs"],
            r["max_abs"],
            q["max_abs"],
            r["mean_sel"],
            q["mean_sel"],
            q["mean_sel"] - r["mean_sel"],
        ])

def overall(rows):
    return {
        "mean_abs": sum(r["abs_residual"] for r in rows)/len(rows),
        "max_abs": max(r["abs_residual"] for r in rows),
    }

def selectivity_roughness(summary_rows):
    ordered = sorted(summary_rows, key=lambda r: r["age"])
    sels = [r["mean_sel"] for r in ordered]
    first_diff_abs_sum = sum(abs(sels[i] - sels[i-1]) for i in range(1, len(sels)))
    second_diff_abs_sum = sum(abs(sels[i] - 2*sels[i-1] + sels[i-2]) for i in range(2, len(sels)))
    return first_diff_abs_sum, second_diff_abs_sum

oa = overall(a)
ob = overall(b)
a_fd, a_sd = selectivity_roughness(sa)
b_fd, b_sd = selectivity_roughness(sb)

lines = []
lines.append("Level 23 smoothness comparison: lambda 0 vs 0.01")
lines.append("=" * 55)
lines.append("")
lines.append("Overall longline residuals")
lines.append("--------------------------")
lines.append(f"lambda 0    mean_abs={oa['mean_abs']:.12g}, max_abs={oa['max_abs']:.12g}")
lines.append(f"lambda 0.01 mean_abs={ob['mean_abs']:.12g}, max_abs={ob['max_abs']:.12g}")
lines.append(f"delta       mean_abs={ob['mean_abs']-oa['mean_abs']:.12g}, max_abs={ob['max_abs']-oa['max_abs']:.12g}")
lines.append("")
lines.append("Longline selectivity roughness")
lines.append("------------------------------")
lines.append(f"lambda 0    first_diff_abs_sum={a_fd:.12g}, second_diff_abs_sum={a_sd:.12g}")
lines.append(f"lambda 0.01 first_diff_abs_sum={b_fd:.12g}, second_diff_abs_sum={b_sd:.12g}")
lines.append(f"delta       first_diff_abs_sum={b_fd-a_fd:.12g}, second_diff_abs_sum={b_sd-a_sd:.12g}")
lines.append("")
lines.append("Decision table")
lines.append("--------------")
lines.append("criterion,value")
lines.append(f"objective_cost_0.01_minus_0,0.1265015")
lines.append(f"mean_abs_residual_cost_0.01_minus_0,{ob['mean_abs']-oa['mean_abs']:.12g}")
lines.append(f"max_abs_residual_cost_0.01_minus_0,{ob['max_abs']-oa['max_abs']:.12g}")
lines.append(f"first_diff_roughness_change_0.01_minus_0,{b_fd-a_fd:.12g}")
lines.append(f"second_diff_roughness_change_0.01_minus_0,{b_sd-a_sd:.12g}")
lines.append("")
lines.append("By-age comparison")
lines.append("-----------------")
lines.append("age,mean_abs_0,mean_abs_0.01,delta_mean_abs,mean_sel_0,mean_sel_0.01,delta_sel")
for r in sa:
    q = b_by_age[r["age"]]
    lines.append(
        f"{r['age']},{r['mean_abs']:.12g},{q['mean_abs']:.12g},"
        f"{q['mean_abs']-r['mean_abs']:.12g},{r['mean_sel']:.12g},"
        f"{q['mean_sel']:.12g},{q['mean_sel']-r['mean_sel']:.12g}"
    )

lines.append("")
lines.append("Interpretation")
lines.append("--------------")
if ob["mean_abs"] <= oa["mean_abs"] * 1.05 and b_sd < a_sd:
    lines.append("lambda=0.01 keeps residual fit essentially unchanged and slightly improves second-difference smoothness.")
elif ob["mean_abs"] <= oa["mean_abs"] * 1.05:
    lines.append("lambda=0.01 keeps residual fit essentially unchanged but does not improve the roughness metric.")
else:
    lines.append("lambda=0.01 worsens residual fit enough to inspect before adopting it.")

txt_path.write_text("\n".join(lines) + "\n")
print(f"wrote: {txt_path}")
print(f"wrote: {csv_path}")
PY

cat "$REPORT"
