#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
import json

IR = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_GRAPH.json")

ir = json.loads(IR.read_text())

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

writers = defaultdict(list)
node_by_key = {}
for node in nodes:
    node_by_key[node["key"]] = node
    for field in node["writes"]:
        writers[field].append(node["key"])

edges = []
seen = set()

def add_edge(source, target, kind, fields=None):
    key = (source, target, kind)
    if source == target or key in seen:
        return
    seen.add(key)
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

for consumer in nodes:
    for field in consumer["reads"]:
        for producer_key in writers.get(field, []):
            producer = node_by_key[producer_key]
            if producer["package"] == consumer["package"]:
                continue
            add_edge(producer_key, consumer["key"], "data_dependency", [field])

def is_external(field):
    return field.startswith("BigeyeModelData") or "Parameters." in field

initial_fields = {"PopulationState.numbers_at_age"}
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
    "graph_version": 1,
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
    for d in diagnostics:
        print(f"  {d['level']}: {d['kind']} {d['operation']} <- {d['field']}")
else:
    print("diagnostics: clean")
