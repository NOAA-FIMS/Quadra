#!/usr/bin/env bash
set -euo pipefail

cat > tools/caa/generate_execution_plan.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path
import json
from collections import defaultdict, deque

IR = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_EXECUTION_PLAN.json")

ir = json.loads(IR.read_text())
packages = ir["packages"]

# Graph scheduler v1
# ------------------
# Builds package dependency edges from state-level consumes/creates/updates.
#
# Important reference-assessment semantics:
# - PopulationState is initialized before AssessmentCycle runs.
# - Packages may create/update the same state.
# - Some packages are allowed to rerun when an input state changes.
#
# The initial PopulationState lets FleetPackage bootstrap FleetState before
# PopulationPackage advances the population. FleetPackage reruns after
# PopulationPackage/MovementPackage updates PopulationState.

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

def state_outputs(pkg):
    return list(dict.fromkeys(pkg["creates"] + pkg["updates"]))

def normalize_phase_name(name):
    return name.replace("Package", "").lower()

def build_static_graph(nodes):
    """Create package dependency graph from state producers to state consumers."""
    producers = defaultdict(list)
    for pkg in nodes:
        for state in state_outputs(pkg):
            producers[state].append(pkg["name"])

    edges = defaultdict(set)
    reasons = defaultdict(list)

    for consumer in nodes:
        cname = consumer["name"]
        for state in state_inputs(consumer):
            for producer in producers.get(state, []):
                if producer == cname:
                    continue
                edges[producer].add(cname)
                reasons[(producer, cname)].append(state)

    return edges, reasons

def topological_layers(nodes, edges):
    names = [p["name"] for p in nodes]
    indegree = {name: 0 for name in names}
    for src, targets in edges.items():
        for dst in targets:
            indegree[dst] += 1

    ready = deque([name for name in names if indegree[name] == 0])
    layers = []

    while ready:
        layer = list(ready)
        ready.clear()
        layers.append(layer)

        for src in layer:
            for dst in sorted(edges.get(src, [])):
                indegree[dst] -= 1
                if indegree[dst] == 0:
                    ready.append(dst)

    unscheduled = [name for name, deg in indegree.items() if deg > 0]
    return layers, unscheduled

pkg_by_name = {p["name"]: p for p in packages}

# For the current reference assessment, the graph has an intentional cycle:
# PopulationPackage consumes FleetState and FleetPackage consumes PopulationState.
# Break that by recognizing externally initialized PopulationState and allowing
# FleetPackage to bootstrap from it before PopulationPackage runs.
bootstrap_nodes = []
main_nodes = []

for pkg in packages:
    if pkg["name"] == "FleetPackage":
        bootstrap_nodes.append(pkg)
    else:
        main_nodes.append(pkg)

phases = []
diagnostics = []
phase_id = 0
available_states = set(INITIAL_STATES)

def append_phase(name, purpose, package_names, derived_from=None):
    global phase_id
    phases.append({
        "phase": len(phases),
        "name": name,
        "purpose": purpose,
        "packages": package_names,
        "derived_from": derived_from or {},
    })

# Phase 0: any package with no state inputs can run before bootstrap.
for pkg in packages:
    if not state_inputs(pkg):
        append_phase(
            normalize_phase_name(pkg["name"]),
            PURPOSES.get(pkg["name"], "execute package"),
            [pkg["name"]],
            {"reason": "no_state_inputs", "consumes": pkg["consumes"], "creates": pkg["creates"], "updates": pkg["updates"]},
        )
        for state in state_outputs(pkg):
            available_states.add(state)

# Bootstrap rerunnable packages whose inputs are already available.
for pkg in bootstrap_nodes:
    missing = [s for s in state_inputs(pkg) if s not in available_states]
    if missing:
        diagnostics.append({
            "level": "error",
            "kind": "bootstrap_missing_inputs",
            "package": pkg["name"],
            "missing": missing,
        })
    else:
        append_phase(
            "fleet_bootstrap" if pkg["name"] == "FleetPackage" else normalize_phase_name(pkg["name"]) + "_bootstrap",
            "initialize fleet mortality before population dynamics" if pkg["name"] == "FleetPackage" else "bootstrap package state",
            [pkg["name"]],
            {"reason": "bootstrap_from_initial_state", "consumes": pkg["consumes"], "creates": pkg["creates"], "updates": pkg["updates"]},
        )
        for state in state_outputs(pkg):
            available_states.add(state)

# Main scheduling: exclude no-input packages already run, then derive graph.
already_run_main = {ph["packages"][0] for ph in phases if ph["name"] != "fleet_bootstrap"}
main_sched_nodes = [
    p for p in packages
    if p["name"] not in already_run_main
]

# Drop the PopulationPackage -> FleetPackage cycle in the main graph by treating
# FleetPackage as rerunnable after PopulationState updates.
edges, reasons = build_static_graph(main_sched_nodes)

# Remove self/semantic bootstrap back-edges where the rerunnable package is the
# producer of an input needed by an upstream package.
for pkg_name in list(edges.keys()):
    for target in list(edges[pkg_name]):
        if pkg_name in RERUNNABLE_PACKAGES and target == "PopulationPackage":
            edges[pkg_name].remove(target)

layers, unscheduled = topological_layers(main_sched_nodes, edges)

if unscheduled:
    diagnostics.append({
        "level": "error",
        "kind": "cycle_or_unscheduled_packages",
        "packages": unscheduled,
    })

for layer in layers:
    # Keep each package in its own phase for now so diagnostics remain clear.
    # Later, packages in the same layer can become a parallel phase.
    for pkg_name in layer:
        pkg = pkg_by_name[pkg_name]
        if pkg_name == "FleetPackage":
            phase_name = "fleet"
            purpose = "recompute fleet predictions after population dynamics"
        else:
            phase_name = normalize_phase_name(pkg_name)
            purpose = PURPOSES.get(pkg_name, "execute package")

        append_phase(
            phase_name,
            purpose,
            [pkg_name],
            {
                "reason": "dependency_graph",
                "consumes": pkg["consumes"],
                "creates": pkg["creates"],
                "updates": pkg["updates"],
                "incoming_edges": [
                    {"from": src, "states": reasons.get((src, pkg_name), [])}
                    for src, targets in edges.items()
                    if pkg_name in targets
                ],
            },
        )

# Validate every package is scheduled at least once.
scheduled = {pkg for ph in phases for pkg in ph["packages"]}
missing_packages = sorted({p["name"] for p in packages} - scheduled)
if missing_packages:
    diagnostics.append({
        "level": "error",
        "kind": "unscheduled_packages",
        "packages": missing_packages,
    })

plan = {
    "name": "Bigeye v2 CAA Execution Plan",
    "source": str(IR),
    "planner_version": 3,
    "scheduler": "graph_dependency_scheduler_v1",
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

echo "updated CAA execution planner to graph-based dependency scheduling"
echo
echo "Run:"
echo "  ./generate_bigeye_v2_caa_execution_plan.sh"
echo "  ./inspect_bigeye_v2_caa_execution_plan.sh"
echo "  ./build_bigeye_v2_caa.sh"
