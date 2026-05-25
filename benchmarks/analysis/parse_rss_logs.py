#!/usr/bin/env python3

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

OUTPUTS = ROOT / "outputs"
NORMALIZED = ROOT / "normalized"

RSS_RE = re.compile(r"Maximum resident set size.*?:\s*([0-9]+)")

def parse_rss_kb(path: Path):
    if not path.exists():
        return ""

    text = path.read_text(errors="replace")
    match = RSS_RE.search(text)

    if not match:
        return ""

    return match.group(1)

def update_csv(csv_path: Path, engine_to_rss: dict):
    if not csv_path.exists():
        return

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    if "peak_rss_kb" not in fieldnames:
        fieldnames.append("peak_rss_kb")

    for row in rows:
        engine = row.get("engine", "")
        rss = engine_to_rss.get(engine, "")

        if rss:
            row["peak_rss_kb"] = rss
        else:
            row.setdefault("peak_rss_kb", "")

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def main():
    rss_map = {
        "quadra": parse_rss_kb(OUTPUTS / "random_intercept_time.log"),
    }

    # Optional future TMB RSS log.
    tmb_rss = parse_rss_kb(OUTPUTS / "tmb_random_intercept_time.log")

    if tmb_rss:
        rss_map["tmb"] = tmb_rss

    update_csv(
        NORMALIZED / "random_intercept_normalized.csv",
        rss_map)

    print("RSS values:")
    for engine, rss in rss_map.items():
        print(f"  {engine}: {rss or 'missing'}")

if __name__ == "__main__":
    main()
