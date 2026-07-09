#!/usr/bin/env bash
set -euo pipefail

cat > tools/caa/generate_ir.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path
import json

ORDER = ["life_history", "population", "movement", "fleet", "observation", "likelihood"]
BASE = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")


def split_list(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def parse_meta(path: Path) -> dict:
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

    out.setdefault("creates", "")
    out.setdefault("updates", "")
    out.setdefault("consumes_fields", "")
    out.setdefault("creates_fields", "")
    out.setdefault("updates_fields", "")
    return out


packages = []
created = {}
updated = {}
consumed = {}
created_fields = {}
updated_fields = {}
consumed_fields = {}

for key in ORDER:
    meta = parse_meta(BASE / key / "package.meta")

    package = {
        "key": key,
        "name": meta["name"],
        "purpose": meta["purpose"],
        "consumes": split_list(meta["consumes"]),
        "produces": split_list(meta["produces"]),
        "creates": split_list(meta["creates"]),
        "updates": split_list(meta["updates"]),
        "consumes_fields": split_list(meta["consumes_fields"]),
        "creates_fields": split_list(meta["creates_fields"]),
        "updates_fields": split_list(meta["updates_fields"]),
        "steps": meta["steps"],
    }

    packages.append(package)

    for item in package["creates"]:
        created.setdefault(item, []).append(package["name"])

    for item in package["updates"]:
        updated.setdefault(item, []).append(package["name"])

    for item in package["consumes"]:
        consumed.setdefault(item, []).append(package["name"])

    for item in package["creates_fields"]:
        created_fields.setdefault(item, []).append(package["name"])

    for item in package["updates_fields"]:
        updated_fields.setdefault(item, []).append(package["name"])

    for item in package["consumes_fields"]:
        consumed_fields.setdefault(item, []).append(package["name"])


diagnostics = []

for state, creators in created.items():
    if len(creators) > 1:
        diagnostics.append({
            "level": "warning",
            "kind": "multiple_creators",
            "state": state,
            "creators": creators,
        })

for field, creators in created_fields.items():
    if len(creators) > 1:
        diagnostics.append({
            "level": "warning",
            "kind": "multiple_field_creators",
            "field": field,
            "creators": creators,
        })

for state, consumers in consumed.items():
    if state.endswith("State") and state not in created and state not in updated:
        diagnostics.append({
            "level": "warning",
            "kind": "missing_creator_or_updater",
            "state": state,
            "consumers": consumers,
        })

for field, consumers in consumed_fields.items():
    if field not in created_fields and field not in updated_fields:
        diagnostics.append({
            "level": "warning",
            "kind": "missing_field_creator_or_updater",
            "field": field,
            "consumers": consumers,
        })

ir = {
    "name": "Bigeye v2 CAA IR",
    "packages": packages,
    "created": created,
    "updated": updated,
    "consumed": consumed,
    "created_fields": created_fields,
    "updated_fields": updated_fields,
    "consumed_fields": consumed_fields,
    "diagnostics": diagnostics,
}

OUT.write_text(json.dumps(ir, indent=2) + "\n")
print(f"wrote {OUT}")

if diagnostics:
    print("diagnostics:")
    for d in diagnostics:
        print(f"  {d['level']}: {d['kind']} {d.get('state', d.get('field', ''))}")
else:
    print("diagnostics: clean")
PY

chmod +x tools/caa/generate_ir.py

python3 - <<'PY'
from pathlib import Path

base = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")

field_meta = {
    "life_history": {
        "consumes_fields": "",
        "creates_fields": "LifeHistoryState.m_at_age, LifeHistoryState.weight_at_age, LifeHistoryState.maturity_at_age",
        "updates_fields": "",
    },
    "population": {
        "consumes_fields": "LifeHistoryState.m_at_age, LifeHistoryState.weight_at_age, LifeHistoryState.maturity_at_age, FleetState.z_at_age",
        "creates_fields": "PopulationState.recruits_by_year, PopulationState.numbers_at_age, PopulationState.survivors_at_age, PopulationState.spawning_biomass_by_year",
        "updates_fields": "PopulationState.recruits_by_year, PopulationState.numbers_at_age, PopulationState.survivors_at_age, PopulationState.spawning_biomass_by_year",
    },
    "movement": {
        "consumes_fields": "PopulationState.numbers_at_age",
        "creates_fields": "",
        "updates_fields": "PopulationState.numbers_at_age",
    },
    "fleet": {
        "consumes_fields": "LifeHistoryState.m_at_age, LifeHistoryState.weight_at_age, PopulationState.numbers_at_age",
        "creates_fields": "FleetState.selectivity_at_age, FleetState.f_at_age, FleetState.z_at_age, FleetState.catch_numbers_at_age, FleetState.catch_biomass_at_age, FleetState.total_catch_biomass_by_year",
        "updates_fields": "FleetState.selectivity_at_age, FleetState.f_at_age, FleetState.z_at_age, FleetState.catch_numbers_at_age, FleetState.catch_biomass_at_age, FleetState.total_catch_biomass_by_year",
    },
    "observation": {
        "consumes_fields": "PopulationState.spawning_biomass_by_year, FleetState.catch_numbers_at_age",
        "creates_fields": "FleetState.predicted_index_by_year, FleetState.predicted_catch_age_proportion",
        "updates_fields": "FleetState.predicted_index_by_year, FleetState.predicted_catch_age_proportion",
    },
    "likelihood": {
        "consumes_fields": "FleetState.total_catch_biomass_by_year, FleetState.predicted_index_by_year, FleetState.predicted_catch_age_proportion",
        "creates_fields": "LikelihoodState.catch_nll, LikelihoodState.index_nll, LikelihoodState.agecomp_nll, LikelihoodState.total_nll",
        "updates_fields": "LikelihoodState.catch_nll, LikelihoodState.index_nll, LikelihoodState.agecomp_nll, LikelihoodState.total_nll",
    },
}

for key, fields in field_meta.items():
    p = base / key / "package.meta"
    lines = p.read_text().splitlines()
    out = []
    for line in lines:
        if line.startswith(("consumes_fields:", "creates_fields:", "updates_fields:")):
            continue
        out.append(line)

    # Insert field metadata after updates.
    final = []
    inserted = False
    for line in out:
        final.append(line)
        if line.startswith("updates:"):
            final.append(f"consumes_fields: {fields['consumes_fields']}")
            final.append(f"creates_fields: {fields['creates_fields']}")
            final.append(f"updates_fields: {fields['updates_fields']}")
            inserted = True
    if not inserted:
        final.append(f"consumes_fields: {fields['consumes_fields']}")
        final.append(f"creates_fields: {fields['creates_fields']}")
        final.append(f"updates_fields: {fields['updates_fields']}")

    p.write_text("\n".join(final) + "\n")
PY

cat > tools/caa/generate_execution_plan.py <<'PY'
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
PY

chmod +x tools/caa/generate_ir.py
chmod +x tools/caa/generate_execution_plan.py

echo "added field-level CAA dependency metadata and scheduler"
echo
echo "Run:"
echo "  ./generate_bigeye_v2_caa_ir.sh"
echo "  ./generate_bigeye_v2_caa_execution_plan.sh"
echo "  ./inspect_bigeye_v2_caa_execution_plan.sh"
echo "  ./build_bigeye_v2_caa.sh"
