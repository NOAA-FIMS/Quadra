#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> int:
    if len(sys.argv) < 3:
        print("usage: make_scaling_plot.py results.csv output.png")
        return 1

    csv_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    rows = []
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    quadra = []
    tmb = []
    quadra_rss = []
    tmb_rss = []

    for row in rows:
        n = int(row["n"])
        q = float(row["quadra_ms"])
        t = float(row["tmb_ms"])
        quadra.append((n, q))
        tmb.append((n, t))
        quadra_rss.append((n, float(row["quadra_peak_rss_mib"])))
        tmb_rss.append((n, float(row["tmb_peak_rss_mib"])))

    quadra.sort()
    tmb.sort()
    quadra_rss.sort()
    tmb_rss.sort()

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot([x for x, _ in quadra], [y for _, y in quadra], marker="o", label="Quadra persistent tridiagonal")
    axes[0].plot([x for x, _ in tmb], [y for _, y in tmb], marker="^", label="TMB AD/Laplace")
    axes[0].set_xlabel("Number of latent-state years")
    axes[0].set_ylabel("Milliseconds per fixed-theta evaluation")
    axes[0].set_yscale("log")
    axes[0].set_title("Runtime")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot([x for x, _ in quadra_rss], [y for _, y in quadra_rss], marker="o", label="Quadra persistent tridiagonal")
    axes[1].plot([x for x, _ in tmb_rss], [y for _, y in tmb_rss], marker="^", label="TMB AD/Laplace")
    axes[1].set_xlabel("Number of latent-state years")
    axes[1].set_ylabel("Peak RSS (MiB)")
    axes[1].set_yscale("log")
    axes[1].set_title("Peak resident memory")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    fig.suptitle("Matched state-space surplus-production Laplace scaling")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)

    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
