#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

def fmt(x):
    x = "" if x is None else str(x)
    if x == "":
        return ""
    try:
        v = float(x)
        if abs(v) >= 100:
            return f"{v:.1f}"
        if abs(v) >= 1:
            return f"{v:.3f}"
        return f"{v:.6f}"
    except ValueError:
        return x

def csv_to_markdown(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        print("_No rows._")
        return

    fields = list(rows[0].keys())

    print("| " + " | ".join(fields) + " |")
    print("| " + " | ".join(["---"] * len(fields)) + " |")

    for row in rows:
        print("| " + " | ".join(fmt(row.get(k, "")) for k in fields) + " |")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: csv_to_markdown.py path/to/file.csv")

    csv_to_markdown(Path(sys.argv[1]))
