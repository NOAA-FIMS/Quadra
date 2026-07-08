#!/usr/bin/env bash
set -euo pipefail

cat > tools/caa/generate_registry.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path

ORDER = ["life_history", "population", "movement", "fleet", "observation", "likelihood"]
BASE = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_PACKAGE_REGISTRY.md")

def parse_meta(path):
    out = {"steps": []}
    in_steps = False
    for line in path.read_text().splitlines():
        line = line.rstrip()
        if line == "steps:":
            in_steps = True
            continue
        if in_steps and line.startswith("  - "):
            out["steps"].append(line[4:])
            continue
        in_steps = False
        if ": " in line:
            k, v = line.split(": ", 1)
            out[k] = v
    return out

rows = [
    "# CAA Package Registry",
    "",
    "Generated from `architecture/packages/*/package.meta`.",
    "",
    "| Package | Consumes | Produces | Steps |",
    "|---|---|---|---|",
]

for d in ORDER:
    meta = parse_meta(BASE / d / "package.meta")
    steps = "<br>".join(meta["steps"])
    rows.append(f"| {meta['name']} | {meta['consumes']} | {meta['produces']} | {steps} |")

OUT.write_text("\n".join(rows) + "\n")
print(f"wrote {OUT}")
PY

chmod +x tools/caa/generate_registry.py

cat > generate_bigeye_v2_caa_registry.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
python3 tools/caa/generate_registry.py
SH

chmod +x generate_bigeye_v2_caa_registry.sh

echo "created CAA registry generator"
