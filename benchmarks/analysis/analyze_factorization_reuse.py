#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

csv_path = ROOT / "exact_laplace_gradient/factorization_reuse_benchmark.csv"
out_path = ROOT / "outputs/factorization_reuse_summary.md"

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

groups = defaultdict(list)

for r in rows:
    groups[r["n_state"]].append(r)

lines = []
lines.append("# Sparse Factorization Reuse Summary")
lines.append("")
lines.append("| n_state | mean_fresh_ms | mean_reuse_ms | mean_reuse_ratio | reuse_improvement_pct | hessian_nnz | factor_nnz | fill_ratio |")
lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")

for n_state in sorted(groups, key=lambda x: int(x)):
    g = groups[n_state]

    fresh = [to_float(r["fresh_ms"]) for r in g if to_float(r["fresh_ms"]) is not None]
    reuse = [to_float(r["reuse_ms"]) for r in g if to_float(r["reuse_ms"]) is not None]
    ratios = [to_float(r["reuse_ratio"]) for r in g if to_float(r["reuse_ratio"]) is not None]

    mean_fresh = sum(fresh) / len(fresh) if fresh else None
    mean_reuse = sum(reuse) / len(reuse) if reuse else None
    mean_ratio = sum(ratios) / len(ratios) if ratios else None
    improvement = 100.0 * (1.0 - mean_ratio) if mean_ratio is not None else None

    hessian_nnz = g[0].get("hessian_nnz", "")
    factor_nnz = g[0].get("factor_nnz", "")
    fill_ratio = g[0].get("fill_ratio", "")

    lines.append(
        f"| {n_state} | {fmt(mean_fresh)} | {fmt(mean_reuse)} | "
        f"{fmt(mean_ratio)} | {fmt(improvement)} | "
        f"{hessian_nnz} | {factor_nnz} | {fill_ratio} |"
    )

lines.append("")
lines.append("Interpretation notes:")
lines.append("")
lines.append("- `fresh_ms` uses a new sparse LDLT object for each factorization.")
lines.append("- `reuse_ms` reuses symbolic analysis and performs numeric refactorization.")
lines.append("- Lower reuse ratios indicate stronger benefit from symbolic reuse.")
lines.append("- This isolates sparse factorization reuse, not full exact-gradient reuse.")

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text("\n".join(lines) + "\n")

print(f"Wrote {out_path}")
