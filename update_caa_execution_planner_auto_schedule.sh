#!/usr/bin/env bash
set -euo pipefail

cat > tools/caa/generate_execution_plan.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path
import json

IR = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_EXECUTION_PLAN.json")

ir = json.loads(IR.read_text())
packages = ir["packages"]

INITIAL_STATES = {"PopulationState"}
RERUNNABLE_PACKAGES = {"FleetPackage"}

PURPOSES = {
    "LifeHistoryPackage": "compute life-history state",
    "PopulationPackage": "advance population state",
    "MovementPackage": "move individuals across populations",
    "FleetPackage": "compute fleet mortality and predictions",
    "ObservationPackage": "predict observations",
    "LikelihoodPackage": "evaluate likelihood components",
}

def state_inputs(pkg):
    return [x for x in pkg["consumes"] if x.endswith("State")]

available = set(INITIAL_STATES)
last_updated = {state: 0 for state in INITIAL_STATES}
last_run = {}
run_count = {}
phases = []
diagnostics = []
phase_id = 0
progress = True

while progress:
    progress = False

    for pkg in packages:
        name = pkg["name"]
        inputs = state_inputs(pkg)
        missing = [s for s in inputs if s not in available]

        if missing:
            continue

        has_run = name in last_run
        stale = any(last_updated.get(s, -1) > last_run.get(name, -1)
                    for s in inputs)

        if has_run:
            if name not in RERUNNABLE_PACKAGES:
                continue
            if not stale:
                continue

        if name in {"ObservationPackage", "LikelihoodPackage"}:
            fleet_runs = run_count.get("FleetPackage", 0)
            population_runs = run_count.get("PopulationPackage", 0)
            if population_runs > 0 and fleet_runs < 2:
                continue

        phase_name = name.replace("Package", "").lower()
        purpose = PURPOSES.get(name, "execute package")

        if name == "FleetPackage" and run_count.get(name, 0) == 0:
            phase_name = "fleet_bootstrap"
            purpose = "initialize fleet mortality before population dynamics"
        elif name == "FleetPackage":
            phase_name = "fleet"
            purpose = "recompute fleet predictions after population dynamics"

        phases.append({
            "phase": phase_id,
            "name": phase_name,
            "purpose": purpose,
            "packages": [name],
            "derived_from": {
                "consumes": pkg["consumes"],
                "creates": pkg["creates"],
                "updates": pkg["updates"],
            },
        })

        phase_id += 1
        progress = True
        last_run[name] = phase_id
        run_count[name] = run_count.get(name, 0) + 1

        for state in pkg["creates"] + pkg["updates"]:
            available.add(state)
            last_updated[state] = phase_id

        break

all_package_names = {p["name"] for p in packages}
planned_once = {pkg for ph in phases for pkg in ph["packages"]}
missing_packages = sorted(all_package_names - planned_once)

if missing_packages:
    diagnostics.append({
        "level": "error",
        "kind": "unscheduled_packages",
        "packages": missing_packages,
    })

for pkg in packages:
    missing = [s for s in state_inputs(pkg) if s not in available]
    if missing:
        diagnostics.append({
            "level": "error",
            "kind": "unsatisfied_state_dependencies",
            "package": pkg["name"],
            "missing": missing,
        })

plan = {
    "name": "Bigeye v2 CAA Execution Plan",
    "source": str(IR),
    "planner_version": 2,
    "scheduler": "state_dependency_scheduler_v1",
    "initial_states": sorted(INITIAL_STATES),
    "rerunnable_packages": sorted(RERUNNABLE_PACKAGES),
    "phases": phases,
    "diagnostics": diagnostics,
}

OUT.write_text(json.dumps(plan, indent=2) + "\n")
print(f"wrote {OUT}")

if diagnostics:
    print("diagnostics:")
    for d in diagnostics:
        print(f"  {d['level']}: {d['kind']}")
    if any(d["level"] == "error" for d in diagnostics):
        raise SystemExit(1)
else:
    print("diagnostics: clean")
PY

chmod +x tools/caa/generate_execution_plan.py

echo "updated CAA execution planner to derive phases from IR state dependencies"
echo
echo "Run:"
echo "  ./generate_bigeye_v2_caa_execution_plan.sh"
echo "  ./inspect_bigeye_v2_caa_execution_plan.sh"
echo "  ./build_bigeye_v2_caa.sh"
