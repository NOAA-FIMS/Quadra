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

    for row in rows:
        n = int(row["n"])
        q = float(row["quadra_ms"])
        t = float(row["tmb_ms"])
        quadra.append((n, q))
        tmb.append((n, t))

    quadra.sort()
    tmb.sort()

    plt.figure()
    plt.plot([x for x, _ in quadra], [y for _, y in quadra], marker="o", label="Quadra analytic tridiagonal")
    plt.plot([x for x, _ in tmb], [y for _, y in tmb], marker="o", label="TMB AD/Laplace")
    plt.xlabel("Number of years / latent states scale")
    plt.ylabel("Milliseconds per fixed-theta Laplace evaluation")
    plt.yscale("log")
    plt.title("State-space surplus production Laplace scaling")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)

    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
