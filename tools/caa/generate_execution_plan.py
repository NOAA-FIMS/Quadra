#!/usr/bin/env python3
from pathlib import Path
import json
from collections import defaultdict, deque

IR = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_EXECUTION_PLAN.json")

ir = json.loads(IR.read_text())
packages = ir["packages"]

INITIAL_FIELDS = {"PopulationState.numbers_at_age"}
RERUNNABLE_PACKAGES = {"FleetPackage"}

PURPOSES = {
    "LifeHistoryPackage": "compute life-history state",
    "PopulationPackage": "advance population state",
    "MovementPackage": "move individuals across populations",
    "FleetPackage": "compute fleet mortality and predictions",
    "ObservationPackage": "predict observations",
    "LikelihoodPackage": "evaluate likelihood components",
}

def inputs(pkg):
    fields = pkg.get("consumes_fields") or []
    if fields:
        return fields
    return [x for x in pkg["consumes"] if x.endswith("State")]

def outputs(pkg):
    fields = (pkg.get("creates_fields") or []) + (pkg.get("updates_fields") or [])
    if fields:
        return list(dict.fromkeys(fields))
    return list(dict.fromkeys(pkg["creates"] + pkg["updates"]))

def normalize_phase_name(name):
    return name.replace("Package", "").lower()

def build_graph(nodes):
    producers = defaultdict(list)
    for pkg in nodes:
        for item in outputs(pkg):
            producers[item].append(pkg["name"])

    edges = defaultdict(set)
    reasons = defaultdict(list)
    for consumer in nodes:
        cname = consumer["name"]
        for item in inputs(consumer):
            for producer in producers.get(item, []):
                if producer == cname:
                    continue
                edges[producer].add(cname)
                reasons[(producer, cname)].append(item)

    return edges, reasons

def topological_layers(nodes, edges):
    names = [p["name"] for p in nodes]
    indegree = {name: 0 for name in names}
    for src, targets in edges.items():
        for dst in targets:
            if dst in indegree:
                indegree[dst] += 1

    ready = deque([name for name in names if indegree[name] == 0])
    layers = []

    while ready:
        layer = list(ready)
        ready.clear()
        layers.append(layer)
        for src in layer:
            for dst in sorted(edges.get(src, [])):
                if dst not in indegree:
                    continue
                indegree[dst] -= 1
                if indegree[dst] == 0:
                    ready.append(dst)

    unscheduled = [name for name, deg in indegree.items() if deg > 0]
    return layers, unscheduled

pkg_by_name = {p["name"]: p for p in packages}
phases = []
diagnostics = []
available = set(INITIAL_FIELDS)

def append_phase(name, purpose, package_names, derived_from=None):
    phases.append({
        "phase": len(phases),
        "name": name,
        "purpose": purpose,
        "packages": package_names,
        "derived_from": derived_from or {},
    })

# Packages with all inputs available can run first, but keep only no-input
# packages here so the bootstrap phase remains explicit.
for pkg in packages:
    if not inputs(pkg):
        append_phase(
            normalize_phase_name(pkg["name"]),
            PURPOSES.get(pkg["name"], "execute package"),
            [pkg["name"]],
            {"reason": "no_inputs", "inputs": inputs(pkg), "outputs": outputs(pkg)},
        )
        available.update(outputs(pkg))

# Bootstrap rerunnable packages whose inputs are already available.
for pkg in packages:
    if pkg["name"] not in RERUNNABLE_PACKAGES:
        continue

    missing = [x for x in inputs(pkg) if x not in available]
    if missing:
        diagnostics.append({
            "level": "error",
            "kind": "bootstrap_missing_inputs",
            "package": pkg["name"],
            "missing": missing,
        })
        continue

    append_phase(
        "fleet_bootstrap" if pkg["name"] == "FleetPackage" else normalize_phase_name(pkg["name"]) + "_bootstrap",
        "initialize fleet mortality before population dynamics" if pkg["name"] == "FleetPackage" else "bootstrap package state",
        [pkg["name"]],
        {"reason": "bootstrap_from_initial_fields", "inputs": inputs(pkg), "outputs": outputs(pkg)},
    )
    available.update(outputs(pkg))

already_run_nonbootstrap = {
    ph["packages"][0]
    for ph in phases
    if not ph["name"].endswith("_bootstrap") and ph["name"] != "fleet_bootstrap"
}

main_nodes = [p for p in packages if p["name"] not in already_run_nonbootstrap]

edges, reasons = build_graph(main_nodes)

# Bootstrap semantics: FleetPackage already produced FleetState.z_at_age before
# PopulationPackage runs, so do not require the main FleetPackage pass to precede
# PopulationPackage.
if "FleetPackage" in edges:
    edges["FleetPackage"].discard("PopulationPackage")

layers, unscheduled = topological_layers(main_nodes, edges)
if unscheduled:
    diagnostics.append({
        "level": "error",
        "kind": "cycle_or_unscheduled_packages",
        "packages": unscheduled,
    })

for layer in layers:
    for pkg_name in layer:
        pkg = pkg_by_name[pkg_name]
        phase_name = normalize_phase_name(pkg_name)
        purpose = PURPOSES.get(pkg_name, "execute package")

        if pkg_name == "FleetPackage":
            phase_name = "fleet"
            purpose = "recompute fleet predictions after population dynamics"

        append_phase(
            phase_name,
            purpose,
            [pkg_name],
            {
                "reason": "field_dependency_graph",
                "inputs": inputs(pkg),
                "outputs": outputs(pkg),
                "incoming_edges": [
                    {"from": src, "items": reasons.get((src, pkg_name), [])}
                    for src, targets in edges.items()
                    if pkg_name in targets
                ],
            },
        )

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
    "planner_version": 4,
    "scheduler": "field_dependency_scheduler_v1",
    "initial_fields": sorted(INITIAL_FIELDS),
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
