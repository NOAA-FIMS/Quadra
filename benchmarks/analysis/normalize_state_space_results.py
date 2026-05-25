#!/usr/bin/env python3

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

quadra_csv = ROOT / "comparisons/tmb_state_space/comparison_outputs/quadra_state_space_compare.csv"
tmb_csv = ROOT / "comparisons/tmb_state_space/comparison_outputs/tmb_state_space_compare.csv"
out_csv = ROOT / "normalized/state_space_normalized.csv"

fieldnames = [
    "engine",
    "model",
    "n_state",
    "n_random",
    "setup_ms",
    "objective_eval_ms",
    "gradient_eval_ms",
    "optimization_ms",
    "workspace_ms",
    "implicit_derivatives_ms",
    "factorization_ms",
    "structure_ms",
    "total_wall_ms",
    "hessian_nnz",
    "hessian_density",
    "factor_nnz",
    "fill_ratio",
    "peak_rss_kb",
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
                "setup_ms": "",
                "objective_eval_ms": r.get("objective_eval_ms", ""),
                "gradient_eval_ms": "",
                "optimization_ms": "",
                "workspace_ms": r.get("workspace_ms", ""),
                "implicit_derivatives_ms": r.get("implicit_derivatives_ms", ""),
                "factorization_ms": r.get("factorization_ms", ""),
                "structure_ms": "",
                "total_wall_ms": r.get("total_wall_ms", ""),
                "hessian_nnz": r.get("hessian_nnz", ""),
                "hessian_density": r.get("hessian_density", ""),
                "factor_nnz": r.get("factor_nnz", ""),
                "fill_ratio": r.get("fill_ratio", ""),
                "peak_rss_kb": "",
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
                "setup_ms": r.get("setup_ms", ""),
                "objective_eval_ms": r.get("objective_eval_ms", ""),
                "gradient_eval_ms": r.get("gradient_eval_ms", ""),
                "optimization_ms": r.get("optimization_ms", ""),
                "workspace_ms": "",
                "implicit_derivatives_ms": "",
                "factorization_ms": "",
                "structure_ms": r.get("structure_ms", ""),
                "total_wall_ms": r.get("total_wall_ms", ""),
                "hessian_nnz": r.get("hessian_nnz", ""),
                "hessian_density": r.get("hessian_density", ""),
                "factor_nnz": r.get("factor_nnz", ""),
                "fill_ratio": r.get("fill_ratio", ""),
                "peak_rss_kb": "",
                "success": r.get("success", "")
            })

out_csv.parent.mkdir(parents=True, exist_ok=True)

with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote normalized state-space benchmark CSV: {out_csv}")
