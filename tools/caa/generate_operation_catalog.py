#!/usr/bin/env python3
from pathlib import Path

ORDER = ["life_history", "population", "movement", "fleet", "observation", "likelihood"]
BASE = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_CATALOG.md")

def parse_meta(path):
    out = {"operations": []}
    in_ops = False
    for line in path.read_text().splitlines():
        line = line.rstrip()
        if line == "operations:":
            in_ops = True
            continue
        if in_ops and line.startswith("  - "):
            out["operations"].append(line[4:])
            continue
        in_ops = False
        if ": " in line:
            k, v = line.split(": ", 1)
            out[k] = v
    return out

lines = [
    "# CAA Operation Catalog",
    "",
    "Generated from `architecture/packages/*/package.meta`.",
    "",
]

for key in ORDER:
    meta = parse_meta(BASE / key / "package.meta")
    lines.extend([
        f"## {meta['name']}",
        "",
        f"**Purpose:** {meta['purpose']}",
        "",
        "**Operations:**",
        "",
    ])

    for op in meta["operations"]:
        lines.append(f"- {op}")

    lines.extend(["", "---", ""])

OUT.write_text("\n".join(lines))
print(f"wrote {OUT}")
