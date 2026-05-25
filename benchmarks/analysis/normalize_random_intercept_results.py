#!/usr/bin/env python3

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

quadra_csv = (
    ROOT /
    "comparisons" /
    "tmb_random_intercept" /
    "comparison_outputs" /
    "quadra_random_intercept_compare.csv"
)

tmb_csv = (
    ROOT /
    "comparisons" /
    "tmb_random_intercept" /
    "comparison_outputs" /
    "tmb_random_intercept_compare.csv"
)

out_csv = (
    ROOT /
    "normalized" /
    "random_intercept_normalized.csv"
)

rows = []

if quadra_csv.exists():
    with open(quadra_csv, newline="") as f:
        reader = csv.DictReader(f)

        for r in reader:
            rows.append({
                "engine": "quadra",
                "model": "random_intercept",
                "n_obs": r.get("n_obs", ""),
                "n_random": 1,
                "workspace_ms": r.get("workspace_ms", ""),
                "implicit_derivatives_ms": r.get("implicit_derivatives_ms", ""),
                "factorization_ms": r.get("factorization_ms", ""),
                "total_wall_ms": r.get("total_wall_ms", ""),
                "peak_rss_kb": "",
                "success": r.get("success", "")
            })

if tmb_csv.exists():
    with open(tmb_csv, newline="") as f:
        reader = csv.DictReader(f)

        for r in reader:
            rows.append({
                "engine": "tmb",
                "model": "random_intercept",
                "n_obs": r.get("n_obs", ""),
                "n_random": 1,
                "workspace_ms": "",
                "implicit_derivatives_ms": "",
                "factorization_ms": "",
                "total_wall_ms": r.get("elapsed_ms", ""),
                "peak_rss_kb": "",
                "success": r.get("convergence", "")
            })

fieldnames = [
    "engine",
    "model",
    "n_obs",
    "n_random",
    "workspace_ms",
    "implicit_derivatives_ms",
    "factorization_ms",
    "total_wall_ms",
    "peak_rss_kb",
    "success"
]

out_csv.parent.mkdir(parents=True, exist_ok=True)

with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()

    for row in rows:
        writer.writerow(row)

print(f"Wrote normalized benchmark CSV: {out_csv}")
