#!/usr/bin/env python3
from pathlib import Path
import json

IR = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_ARCHITECTURE_DIAGRAM.md")

ir = json.loads(IR.read_text())
packages = ir["packages"]

lines = [
    "# CAA Architecture Diagram",
    "",
    "Generated from `CAA_IR.json`.",
    "",
    "```mermaid",
    "flowchart TD",
    "  AssessmentCycle[AssessmentCycle]",
]

previous = "AssessmentCycle"

for pkg in packages:
    name = pkg["name"]
    lines.append(f"  {previous} --> {name}")
    previous = name

    for step in pkg["steps"]:
        lines.append(f"  {name} --> {step}")

lines.extend([
    "```",
    "",
])

OUT.write_text("\n".join(lines))
print(f"wrote {OUT}")
