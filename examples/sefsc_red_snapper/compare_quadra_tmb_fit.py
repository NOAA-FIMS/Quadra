#!/usr/bin/env python3
from pathlib import Path
import csv
import math

out = Path("examples/sefsc_red_snapper/outputs")

def read_summary(path):
    d = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                d[row["field"]] = float(row["value"])
            except Exception:
                d[row["field"]] = row["value"]
    return d

q = read_summary(out / "quadra_fit_summary.csv")
t = read_summary(out / "tmb_fit_summary.csv")

fields = ["objective", "r0", "fbar", "q", "sel_a50", "sel_slope", "random_effects"]
path = out / "quadra_vs_tmb_fit_comparison.csv"

with open(path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["field", "quadra", "tmb", "difference", "relative_difference"])
    for field in fields:
        qv = q.get(field, "")
        tv = t.get(field, "")
        diff = ""
        rel = ""
        if isinstance(qv, float) and isinstance(tv, float):
            diff = qv - tv
            rel = diff / tv if tv != 0 and math.isfinite(tv) else ""
        w.writerow([field, qv, tv, diff, rel])

print(f"wrote: {path}")
