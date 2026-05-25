#!/usr/bin/env python3

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

quadra_csv = ROOT / "exact_laplace_gradient/state_space_exact_gradient_benchmark.csv"
tmb_csv = ROOT / "comparisons/tmb_state_space/comparison_outputs/tmb_state_space_compare.csv"
out_csv = ROOT / "normalized/exact_gradient_state_space_comparison.csv"

fieldnames = [
    "engine",
    "model",
    "n_state",
    "n_random",
    "exact_gradient_ms",
    "objective_eval_ms",
    "gradient_eval_ms",
    "fn_gr_eval_ms",
    "hessian_nnz",
    "factor_nnz",
    "fill_ratio",
    "success"
]

rows = []

if quadra_csv.exists():
    with open(quadra_csv, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "engine": "quadra",
                "model": "state_space",
                "n_state": r.get("n_state", ""),
                "n_random": r.get("n_random", ""),
                "exact_gradient_ms": r.get("total_ms", ""),
                "objective_eval_ms": "",
                "gradient_eval_ms": "",
                "fn_gr_eval_ms": "",
                "hessian_nnz": r.get("hessian_nnz", ""),
                "factor_nnz": r.get("factor_nnz", ""),
                "fill_ratio": r.get("fill_ratio", ""),
                "success": r.get("success", "")
            })

if tmb_csv.exists():
    with open(tmb_csv, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "engine": "tmb",
                "model": "state_space",
                "n_state": r.get("n_state", ""),
                "n_random": r.get("n_random", ""),
                "exact_gradient_ms": "",
                "objective_eval_ms": r.get("objective_eval_ms", ""),
                "gradient_eval_ms": r.get("gradient_eval_ms", ""),
                "fn_gr_eval_ms": r.get("fn_gr_eval_ms", ""),
                "hessian_nnz": r.get("hessian_nnz", ""),
                "factor_nnz": r.get("factor_nnz", ""),
                "fill_ratio": r.get("fill_ratio", ""),
                "success": r.get("success", "")
            })

out_csv.parent.mkdir(parents=True, exist_ok=True)

with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {out_csv}")
