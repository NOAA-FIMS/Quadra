#!/usr/bin/env bash
set -euo pipefail

cat > tools/caa/validate_operation_graph.py <<'PY'
#!/usr/bin/env python3
from collections import Counter, defaultdict
from pathlib import Path
import json
import sys

GRAPH = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_GRAPH.json")
PLAN = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_EXECUTION_PLAN.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_GRAPH_VALIDATION.json")

graph = json.loads(GRAPH.read_text())
plan = json.loads(PLAN.read_text()) if PLAN.exists() else {}

nodes = graph["nodes"]
edges = graph["edges"]
initial_fields = set(graph.get("initial_fields", []))
rerunnable_packages = set(plan.get("rerunnable_packages", []))

errors = []
warnings = []


def add_error(kind, **details):
    errors.append({"level": "error", "kind": kind, **details})


def add_warning(kind, **details):
    warnings.append({"level": "warning", "kind": kind, **details})


# ---------------------------------------------------------------------------
# Identity and schema validation
# ---------------------------------------------------------------------------
keys = [node["key"] for node in nodes]
node_keys = set(keys)
node_by_key = {node["key"]: node for node in nodes}

for key, count in Counter(keys).items():
    if count > 1:
        add_error("duplicate_operation_key", operation=key, count=count)

ids_by_package = defaultdict(list)
orders_by_package = defaultdict(list)

for node in nodes:
    for required in ("key", "id", "package", "order", "reads", "writes"):
        if required not in node:
            add_error(
                "missing_operation_attribute",
                operation=node.get("key", "(unknown)"),
                attribute=required,
            )

    package = node.get("package", "(unknown)")
    ids_by_package[package].append(node.get("id"))
    orders_by_package[package].append((node.get("order"), node.get("key")))

    if not isinstance(node.get("order"), int) or node.get("order", -1) < 0:
        add_error(
            "invalid_operation_order",
            operation=node.get("key", "(unknown)"),
            order=node.get("order"),
        )

    if not isinstance(node.get("reads"), list):
        add_error("invalid_reads_type", operation=node.get("key", "(unknown)"))

    if not isinstance(node.get("writes"), list):
        add_error("invalid_writes_type", operation=node.get("key", "(unknown)"))

for package, ids in ids_by_package.items():
    for operation_id, count in Counter(ids).items():
        if count > 1:
            add_error(
                "duplicate_operation_id",
                package=package,
                operation=operation_id,
                count=count,
            )

for package, ordered in orders_by_package.items():
    counts = Counter(order for order, _ in ordered)
    for order, count in counts.items():
        if count > 1:
            add_error(
                "duplicate_declared_order",
                package=package,
                order=order,
                operations=[key for value, key in ordered if value == order],
            )


# ---------------------------------------------------------------------------
# Field lineage validation
# ---------------------------------------------------------------------------
writers = defaultdict(list)
readers = defaultdict(list)

for node in nodes:
    for field in node.get("writes", []):
        writers[field].append(node["key"])
    for field in node.get("reads", []):
        readers[field].append(node["key"])


def is_external(field):
    return (
        field.startswith("BigeyeModelData")
        or "Parameters." in field
        or field.endswith("Parameters")
    )


for field, field_writers in sorted(writers.items()):
    unique_writers = sorted(set(field_writers))
    if len(unique_writers) > 1:
        add_warning(
            "multiple_operation_writers",
            field=field,
            writers=unique_writers,
        )

for field, field_readers in sorted(readers.items()):
    if field in writers or field in initial_fields or is_external(field):
        continue
    add_error(
        "operation_read_without_writer",
        field=field,
        readers=sorted(set(field_readers)),
    )

for field, field_writers in sorted(writers.items()):
    if field in readers or field == "LikelihoodState.total_nll":
        continue
    add_warning(
        "operation_output_never_consumed",
        field=field,
        writers=sorted(set(field_writers)),
    )


# ---------------------------------------------------------------------------
# Edge validation
# ---------------------------------------------------------------------------
seen_edges = set()
adjacency = defaultdict(set)
edge_by_pair = defaultdict(list)

for edge in edges:
    source = edge.get("from")
    target = edge.get("to")
    kind = edge.get("kind")

    if source not in node_keys:
        add_error("edge_source_missing", edge=edge)
    if target not in node_keys:
        add_error("edge_target_missing", edge=edge)
    if source == target:
        add_error("operation_self_edge", edge=edge)
    if kind not in {"declared_order", "data_dependency"}:
        add_error("unknown_operation_edge_kind", edge=edge)

    signature = (source, target, kind)
    if signature in seen_edges:
        add_warning("duplicate_operation_edge", edge=edge)
    seen_edges.add(signature)

    if source in node_keys and target in node_keys:
        adjacency[source].add(target)
        edge_by_pair[(source, target)].append(edge)


# ---------------------------------------------------------------------------
# Strongly connected components
# ---------------------------------------------------------------------------
index = 0
stack = []
on_stack = set()
indices = {}
lowlink = {}
components = []


def strongconnect(node):
    global index

    indices[node] = index
    lowlink[node] = index
    index += 1
    stack.append(node)
    on_stack.add(node)

    for target in sorted(adjacency.get(node, [])):
        if target not in indices:
            strongconnect(target)
            lowlink[node] = min(lowlink[node], lowlink[target])
        elif target in on_stack:
            lowlink[node] = min(lowlink[node], indices[target])

    if lowlink[node] == indices[node]:
        component = []
        while True:
            member = stack.pop()
            on_stack.remove(member)
            component.append(member)
            if member == node:
                break
        components.append(sorted(component))


for node in sorted(node_keys):
    if node not in indices:
        strongconnect(node)

for component in components:
    has_self_loop = (
        len(component) == 1
        and component[0] in adjacency.get(component[0], set())
    )

    if len(component) == 1 and not has_self_loop:
        continue

    component_set = set(component)
    packages = sorted({
        node_by_key[key]["package"]
        for key in component
        if key in node_by_key
    })

    feedback_fields = set()
    for source in component:
        for target in adjacency.get(source, []):
            if target not in component_set:
                continue
            for edge in edge_by_pair[(source, target)]:
                feedback_fields.update(edge.get("fields", []))

    resolving_packages = sorted(set(packages) & rerunnable_packages)

    details = {
        "operations": component,
        "packages": packages,
        "fields": sorted(feedback_fields),
        "rerunnable_packages": resolving_packages,
        "resolution": (
            "execution planner bootstrap/rerun"
            if resolving_packages
            else "unresolved"
        ),
    }

    if resolving_packages:
        add_warning("operation_graph_feedback_cycle", **details)
    else:
        add_error("operation_graph_cycle", **details)


# ---------------------------------------------------------------------------
# Declared-order contradiction validation
# ---------------------------------------------------------------------------
for edge in edges:
    if edge.get("kind") != "data_dependency":
        continue

    source = node_by_key.get(edge.get("from"))
    target = node_by_key.get(edge.get("to"))
    if source is None or target is None:
        continue

    if (
        source["package"] == target["package"]
        and source["order"] > target["order"]
    ):
        add_error(
            "data_dependency_contradicts_declared_order",
            source=source["key"],
            source_order=source["order"],
            target=target["key"],
            target_order=target["order"],
            fields=edge.get("fields", []),
        )


result = {
    "name": "Bigeye v2 CAA Operation Graph Validation",
    "source": str(GRAPH),
    "execution_plan": str(PLAN) if PLAN.exists() else None,
    "errors": errors,
    "warnings": warnings,
}

OUT.write_text(json.dumps(result, indent=2) + "\n")

print("CAA Operation Graph Validation")
print(f"nodes: {len(nodes)}")
print(f"edges: {len(edges)}")
print(f"errors: {len(errors)}")
print(f"warnings: {len(warnings)}")
print()

for diagnostic in errors + warnings:
    print(f"{diagnostic['level'].upper()}: {diagnostic['kind']}")
    for key, value in diagnostic.items():
        if key not in {"level", "kind"}:
            print(f"  {key}: {value}")
    print()

print(f"wrote {OUT}")

if errors:
    print("FAILED: CAA operation graph validation")
    sys.exit(1)

print("PASSED: CAA operation graph validation")
PY

chmod +x tools/caa/validate_operation_graph.py

cat > validate_bigeye_v2_caa_operation_graph.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

./generate_bigeye_v2_caa_operation_graph.sh
./generate_bigeye_v2_caa_execution_plan.sh
python3 tools/caa/validate_operation_graph.py
SH

chmod +x validate_bigeye_v2_caa_operation_graph.sh

echo "replaced operation-graph validator with planner-aware SCC validation"
echo
echo "Run:"
echo "  ./validate_bigeye_v2_caa_operation_graph.sh"
echo "  ./build_bigeye_v2_caa.sh"
echo "  ./run_bigeye_v2_regression_suite.sh"
