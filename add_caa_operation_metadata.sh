#!/usr/bin/env bash
set -euo pipefail

echo "Adding CAA operation-level metadata..."

BASE="examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages"

python3 - <<'PY'
from pathlib import Path

base = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")

ops = {
    "life_history": [
        "ComputeNaturalMortality",
        "ComputeWeightAtAge",
        "ComputeMaturityAtAge",
    ],
    "population": [
        "Recruitment",
        "Survival",
        "Aging",
        "PlusGroup",
        "SpawningBiomass",
    ],
    "movement": [
        "IdentityMovement",
    ],
    "fleet": [
        "LogisticSelectivity",
        "FishingMortality",
        "BaranovCatch",
    ],
    "observation": [
        "BiomassIndexPrediction",
        "CatchAgeCompositionPrediction",
    ],
    "likelihood": [
        "LognormalCatchLikelihood",
        "LognormalIndexLikelihood",
        "MultinomialAgeCompLikelihood",
    ],
}

for key, operations in ops.items():
    p = base / key / "package.meta"
    lines = p.read_text().splitlines()

    out = []
    in_ops = False
    for line in lines:
        if line == "operations:":
            in_ops = True
            continue
        if in_ops and line.startswith("  - "):
            continue
        in_ops = False
        out.append(line)

    out.append("operations:")
    for op in operations:
        out.append(f"  - {op}")

    p.write_text("\n".join(out) + "\n")
PY

cat > tools/caa/generate_operation_catalog.py <<'PY'
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
PY

chmod +x tools/caa/generate_operation_catalog.py

cat > generate_bigeye_v2_caa_operation_catalog.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

python3 tools/caa/generate_operation_catalog.py
SH

chmod +x generate_bigeye_v2_caa_operation_catalog.sh

python3 - <<'PY'
from pathlib import Path

p = Path("tools/caa/generate_ir.py")
s = p.read_text()

if '"operations": meta.get("operations", [])' not in s:
    s = s.replace(
        'out = {"steps": []}',
        'out = {"steps": [], "operations": []}'
    )

    s = s.replace(
        'in_steps = False',
        'in_steps = False\n    in_operations = False',
        1
    )

    s = s.replace(
'''        if line == "steps:":
            in_steps = True
            continue

        if in_steps and line.startswith("  - "):
            out["steps"].append(line[4:])
            continue

        in_steps = False''',
'''        if line == "steps:":
            in_steps = True
            in_operations = False
            continue

        if line == "operations:":
            in_operations = True
            in_steps = False
            continue

        if in_steps and line.startswith("  - "):
            out["steps"].append(line[4:])
            continue

        if in_operations and line.startswith("  - "):
            out["operations"].append(line[4:])
            continue

        in_steps = False
        in_operations = False'''
    )

    s = s.replace(
'''        "steps": meta["steps"],''',
'''        "steps": meta["steps"],
        "operations": meta.get("operations", []),'''
    )

p.write_text(s)
PY

python3 - <<'PY'
from pathlib import Path

p = Path("build_bigeye_v2_caa.sh")
s = p.read_text()

if "== Generate operation catalog ==" not in s:
    s = s.replace(
'''echo
echo "== Generate registry =="
./generate_bigeye_v2_caa_registry.sh''',
'''echo
echo "== Generate registry =="
./generate_bigeye_v2_caa_registry.sh

echo
echo "== Generate operation catalog =="
./generate_bigeye_v2_caa_operation_catalog.sh'''
    )

p.write_text(s)
PY

./generate_bigeye_v2_caa_operation_catalog.sh
./generate_bigeye_v2_caa_ir.sh

echo
echo "Added CAA operation-level metadata."
echo
echo "Run:"
echo "  head -80 examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_CATALOG.md"
echo "  ./build_bigeye_v2_caa.sh"
echo "  ./run_bigeye_v2_regression_suite.sh"
