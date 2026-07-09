#!/usr/bin/env python3
from pathlib import Path
import json

IR = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_EXECUTION_PLAN.json")

ir = json.loads(IR.read_text())
available = {p["name"] for p in ir["packages"]}

candidate_phases = [
    {"phase": 0, "name": "life_history", "purpose": "compute life-history state", "packages": ["LifeHistoryPackage"]},
    {"phase": 1, "name": "fleet_bootstrap", "purpose": "initialize fleet mortality before population dynamics", "packages": ["FleetPackage"]},
    {"phase": 2, "name": "population", "purpose": "advance population state", "packages": ["PopulationPackage"]},
    {"phase": 3, "name": "movement", "purpose": "move individuals across populations", "packages": ["MovementPackage"]},
    {"phase": 4, "name": "fleet", "purpose": "recompute fleet predictions after population dynamics", "packages": ["FleetPackage"]},
    {"phase": 5, "name": "observation", "purpose": "predict observations", "packages": ["ObservationPackage"]},
    {"phase": 6, "name": "likelihood", "purpose": "evaluate likelihood components", "packages": ["LikelihoodPackage"]},
]

phases = []
diagnostics = []

for phase in candidate_phases:
    missing = [pkg for pkg in phase["packages"] if pkg not in available]
    if missing:
        diagnostics.append({
            "level": "warning",
            "kind": "missing_package_for_phase",
            "phase": phase["name"],
            "missing": missing,
        })
        continue
    phases.append(phase)

plan = {
    "name": "Bigeye v2 CAA Execution Plan",
    "source": str(IR),
    "planner_version": 1,
    "phases": phases,
    "diagnostics": diagnostics,
}

OUT.write_text(json.dumps(plan, indent=2) + "\n")
print(f"wrote {OUT}")

if diagnostics:
    print("diagnostics:")
    for d in diagnostics:
        print(f"  {d['level']}: {d['kind']} {d.get('phase', '')}")
else:
    print("diagnostics: clean")
