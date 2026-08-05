#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path


QUADRA_ROW = re.compile(
    r"^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
    r"\s+(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
    r"\s+([-+0-9.eE]+)\s*$"
)
TMB_ROW = re.compile(r"^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$")
RSS_ROW = re.compile(r"^PEAK_RSS_BYTES n=(\d+) value=(\d+)$")


def parse(path: Path, row_pattern: re.Pattern[str], is_quadra: bool):
    rows = {}
    rss = {}
    for line in path.read_text().splitlines():
        rss_match = RSS_ROW.match(line)
        if rss_match:
            rss[int(rss_match.group(1))] = int(rss_match.group(2))
            continue
        match = row_pattern.match(line)
        if not match:
            continue
        n = int(match.group(1))
        if is_quadra:
            rows[n] = {
                "objective": float(match.group(2)),
                "joint": float(match.group(3)),
                "logdet": float(match.group(4)),
                "nnz": int(match.group(5)),
                "grad_norm": float(match.group(6)),
                "analytic_ms": float(match.group(7)),
                "ms": float(match.group(8)),
            }
        else:
            rows[n] = {"objective": float(match.group(2)), "ms": float(match.group(3))}
    for n, value in rss.items():
        if n in rows:
            rows[n]["peak_rss_mib"] = value / (1024.0 * 1024.0)
    return rows


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: normalize_results.py quadra.log tmb.log results.csv")
        return 1

    quadra = parse(Path(sys.argv[1]), QUADRA_ROW, True)
    tmb = parse(Path(sys.argv[2]), TMB_ROW, False)
    ns = sorted(set(quadra) & set(tmb))
    if not ns:
        raise RuntimeError("no matching Quadra/TMB benchmark rows")

    fields = [
        "n", "quadra_objective", "tmb_objective", "objective_diff",
        "quadra_ms", "tmb_ms", "quadra_speedup",
        "quadra_peak_rss_mib", "tmb_peak_rss_mib", "rss_ratio",
        "joint", "logdet", "nnz", "grad_norm", "quadra_analytic_ms",
    ]
    with Path(sys.argv[3]).open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for n in ns:
            q, t = quadra[n], tmb[n]
            q_rss, t_rss = q["peak_rss_mib"], t["peak_rss_mib"]
            writer.writerow({
                "n": n,
                "quadra_objective": q["objective"],
                "tmb_objective": t["objective"],
                "objective_diff": q["objective"] - t["objective"],
                "quadra_ms": q["ms"],
                "tmb_ms": t["ms"],
                "quadra_speedup": t["ms"] / q["ms"],
                "quadra_peak_rss_mib": q_rss,
                "tmb_peak_rss_mib": t_rss,
                "rss_ratio": t_rss / q_rss,
                "joint": q["joint"],
                "logdet": q["logdet"],
                "nnz": q["nnz"],
                "grad_norm": q["grad_norm"],
                "quadra_analytic_ms": q["analytic_ms"],
            })
    print(f"wrote {sys.argv[3]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
