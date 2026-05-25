#!/usr/bin/env python3

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

csv_path = ROOT / "exact_laplace_gradient/exact_gradient_reuse_benchmark.csv"
out_path = ROOT / "outputs/exact_gradient_reuse_summary.md"

def to_float(x):
    try:
        return float(x)
    except Exception:
        return None

def fmt(x):
    if x is None:
        return "NA"
    if abs(x) >= 100:
        return f"{x:.1f}"
    if abs(x) >= 1:
        return f"{x:.3f}"
    return f"{x:.6f}"

with open(csv_path, newline="") as f:
    rows = list(csv.DictReader(f))

if not rows:
    raise SystemExit("No reuse benchmark rows found.")

metrics = [
    "total_ms",
    "objective_ms",
    "tape_setup_ms",
    "reverse_pass_ms",
    "gradient_extract_ms",
    "total_gradient_ms",
]

cold = rows[0]
warm_rows = rows[max(1, len(rows) // 2):]

lines = []
lines.append("# Exact Laplace Gradient Reuse Summary")
lines.append("")
lines.append("This summary compares the first evaluation against the mean of later evaluations.")
lines.append("")
lines.append("| metric | cold_ms | warm_mean_ms | warm/cold | improvement_pct |")
lines.append("|---|---:|---:|---:|---:|")

for metric in metrics:
    cold_value = to_float(cold.get(metric))
    warm_values = [
        to_float(r.get(metric))
        for r in warm_rows
        if to_float(r.get(metric)) is not None
    ]

    if cold_value is None or cold_value == 0.0 or not warm_values:
        warm_mean = None
        ratio = None
        improvement = None
    else:
        warm_mean = sum(warm_values) / len(warm_values)
        ratio = warm_mean / cold_value
        improvement = 100.0 * (1.0 - ratio)

    lines.append(
        f"| {metric} | {fmt(cold_value)} | {fmt(warm_mean)} | "
        f"{fmt(ratio)} | {fmt(improvement)} |"
    )

lines.append("")
lines.append("Interpretation notes:")
lines.append("")
lines.append("- Values below 1.0 in the warm/cold column indicate faster later evaluations.")
lines.append("- Positive improvement percentages indicate apparent warmup or reuse benefit.")
lines.append("- This benchmark perturbs theta slightly between iterations, so it is a repeated-evaluation workload rather than a literal identical-call cache test.")

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text("\n".join(lines) + "\n")

print(f"Wrote {out_path}")
