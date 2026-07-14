#!/usr/bin/env bash
set -euo pipefail

cat > tools/caa/generate_operation_graph.py <<'PY'
#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
import json

IR = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")
PLAN = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_EXECUTION_PLAN.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_GRAPH.json")

ir = json.loads(IR.read_text())
plan = json.loads(PLAN.read_text()) if PLAN.exists() else {}

nodes = []
for package in ir["packages"]:
    for operation in package.get("operation_nodes", []):
        nodes.append({
            "key": f"{package['name']}::{operation['id']}",
            "id": operation["id"],
            "package": package["name"],
            "order": operation["order"],
            "reads": operation["reads"],
            "writes": operation["writes"],
        })

node_by_key = {node["key"]: node for node in nodes}

package_phases = defaultdict(list)
for phase in plan.get("phases", []):
    for package in phase.get("packages", []):
        package_phases[package].append(phase["phase"])

writers = defaultdict(list)
for node in nodes:
    for field in node["writes"]:
        writers[field].append(node["key"])

edges = []
seen = set()

def add_edge(source, target, kind, fields=None):
    signature = (source, target, kind)
    if source == target or signature in seen:
        return
    seen.add(signature)
    edges.append({
        "from": source,
        "to": target,
        "kind": kind,
        "fields": sorted(set(fields or [])),
    })

by_package = defaultdict(list)
for node in nodes:
    by_package[node["package"]].append(node)

for package_nodes in by_package.values():
    ordered = sorted(package_nodes, key=lambda node: (node["order"], node["id"]))
    for left, right in zip(ordered, ordered[1:]):
        add_edge(left["key"], right["key"], "declared_order")

def can_precede(producer_package, consumer_package):
    producer_positions = package_phases.get(producer_package, [])
    consumer_positions = package_phases.get(consumer_package, [])

    if not producer_positions or not consumer_positions:
        return True

    return any(
        producer_phase < consumer_phase
        for producer_phase in producer_positions
        for consumer_phase in consumer_positions
    )

for consumer in nodes:
    for field in consumer["reads"]:
        for producer_key in writers.get(field, []):
            producer = node_by_key[producer_key]

            if producer["package"] == consumer["package"]:
                continue

            if not can_precede(producer["package"], consumer["package"]):
                continue

            add_edge(
                producer_key,
                consumer["key"],
                "data_dependency",
                [field],
            )

def is_external(field):
    return (
        field.startswith("BigeyeModelData")
        or "Parameters." in field
        or field.endswith("Parameters")
    )

initial_fields = set(plan.get("initial_fields", []))
initial_fields.add("PopulationState.numbers_at_age")

diagnostics = []

for node in nodes:
    for field in node["reads"]:
        if field in writers or field in initial_fields or is_external(field):
            continue
        diagnostics.append({
            "level": "warning",
            "kind": "operation_read_without_writer",
            "operation": node["key"],
            "field": field,
        })

graph = {
    "name": "Bigeye v2 CAA Operation Graph",
    "source": str(IR),
    "execution_plan": str(PLAN) if PLAN.exists() else None,
    "graph_version": 2,
    "dependency_model": "plan_directed_field_lineage_v1",
    "initial_fields": sorted(initial_fields),
    "nodes": nodes,
    "edges": edges,
    "diagnostics": diagnostics,
}

OUT.write_text(json.dumps(graph, indent=2) + "\n")
print(f"wrote {OUT}")
print(f"nodes: {len(nodes)}")
print(f"edges: {len(edges)}")

if diagnostics:
    print("diagnostics:")
    for diagnostic in diagnostics:
        print(
            f"  {diagnostic['level']}: {diagnostic['kind']} "
            f"{diagnostic['operation']} <- {diagnostic['field']}"
        )
else:
    print("diagnostics: clean")
PY

chmod +x tools/caa/generate_operation_graph.py

cat > generate_bigeye_v2_caa_operation_graph.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
./generate_bigeye_v2_caa_ir.sh
./generate_bigeye_v2_caa_execution_plan.sh
python3 tools/caa/generate_operation_graph.py
SH

chmod +x generate_bigeye_v2_caa_operation_graph.sh

echo "updated operation graph to use execution-plan-directed field lineage"
echo
echo "Run:"
echo "  ./validate_bigeye_v2_caa_operation_graph.sh"
echo "  ./build_bigeye_v2_caa.sh"
echo "  ./run_bigeye_v2_regression_suite.sh"
