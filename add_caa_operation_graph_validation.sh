#!/usr/bin/env bash
set -euo pipefail

cat > tools/caa/validate_operation_graph.py <<'PY'
#!/usr/bin/env python3
from collections import Counter, defaultdict, deque
from pathlib import Path
import json
import sys

GRAPH = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_GRAPH.json")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_GRAPH_VALIDATION.json")

graph = json.loads(GRAPH.read_text())
nodes = graph["nodes"]
edges = graph["edges"]
initial_fields = set(graph.get("initial_fields", []))

errors = []
warnings = []

def error(kind, **details):
    errors.append({"level": "error", "kind": kind, **details})

def warning(kind, **details):
    warnings.append({"level": "warning", "kind": kind, **details})

keys = [node["key"] for node in nodes]
node_keys = set(keys)
node_by_key = {node["key"]: node for node in nodes}

for key, count in Counter(keys).items():
    if count > 1:
        error("duplicate_operation_key", operation=key, count=count)

ids_by_package = defaultdict(list)
orders_by_package = defaultdict(list)

for node in nodes:
    for required in ("key", "id", "package", "order", "reads", "writes"):
        if required not in node:
            error("missing_operation_attribute",
                  operation=node.get("key", "(unknown)"),
                  attribute=required)

    ids_by_package[node.get("package", "(unknown)")].append(node.get("id"))
    orders_by_package[node.get("package", "(unknown)")].append(
        (node.get("order"), node.get("key"))
    )

    if not isinstance(node.get("order"), int) or node.get("order", -1) < 0:
        error("invalid_operation_order",
              operation=node.get("key", "(unknown)"),
              order=node.get("order"))

    if not isinstance(node.get("reads"), list):
        error("invalid_reads_type", operation=node.get("key", "(unknown)"))

    if not isinstance(node.get("writes"), list):
        error("invalid_writes_type", operation=node.get("key", "(unknown)"))

for package, ids in ids_by_package.items():
    for operation_id, count in Counter(ids).items():
        if count > 1:
            error("duplicate_operation_id",
                  package=package,
                  operation=operation_id,
                  count=count)

for package, ordered in orders_by_package.items():
    counts = Counter(order for order, _ in ordered)
    for order, count in counts.items():
        if count > 1:
            error("duplicate_declared_order",
                  package=package,
                  order=order,
                  operations=[key for value, key in ordered if value == order])

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
    unique = sorted(set(field_writers))
    if len(unique) > 1:
        warning("multiple_operation_writers", field=field, writers=unique)

for field, field_readers in sorted(readers.items()):
    if field in writers or field in initial_fields or is_external(field):
        continue
    error("operation_read_without_writer",
          field=field,
          readers=sorted(set(field_readers)))

for field, field_writers in sorted(writers.items()):
    if field in readers or field == "LikelihoodState.total_nll":
        continue
    warning("operation_output_never_consumed",
            field=field,
            writers=sorted(set(field_writers)))

seen_edges = set()
adjacency = defaultdict(set)
indegree = {key: 0 for key in node_keys}

for edge in edges:
    source = edge.get("from")
    target = edge.get("to")
    kind = edge.get("kind")

    if source not in node_keys:
        error("edge_source_missing", edge=edge)
    if target not in node_keys:
        error("edge_target_missing", edge=edge)
    if source == target:
        error("operation_self_edge", edge=edge)
    if kind not in {"declared_order", "data_dependency"}:
        error("unknown_operation_edge_kind", edge=edge)

    signature = (source, target, kind)
    if signature in seen_edges:
        warning("duplicate_operation_edge", edge=edge)
    seen_edges.add(signature)

    if source in node_keys and target in node_keys and target not in adjacency[source]:
        adjacency[source].add(target)
        indegree[target] += 1

ready = deque(sorted(key for key, degree in indegree.items() if degree == 0))
visited = []

while ready:
    key = ready.popleft()
    visited.append(key)
    for target in sorted(adjacency.get(key, [])):
        indegree[target] -= 1
        if indegree[target] == 0:
            ready.append(target)

if len(visited) != len(node_keys):
    error("operation_graph_cycle",
          operations=sorted(key for key, degree in indegree.items()
                            if degree > 0))

for edge in edges:
    if edge.get("kind") != "data_dependency":
        continue

    source = node_by_key.get(edge.get("from"))
    target = node_by_key.get(edge.get("to"))
    if source is None or target is None:
        continue

    if source["package"] == target["package"] and source["order"] > target["order"]:
        error("data_dependency_contradicts_declared_order",
              source=source["key"],
              source_order=source["order"],
              target=target["key"],
              target_order=target["order"],
              fields=edge.get("fields", []))

result = {
    "name": "Bigeye v2 CAA Operation Graph Validation",
    "source": str(GRAPH),
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
python3 tools/caa/validate_operation_graph.py
SH

chmod +x validate_bigeye_v2_caa_operation_graph.sh

python3 - <<'PY'
from pathlib import Path

p = Path("build_bigeye_v2_caa.sh")
s = p.read_text()

if "== Validate operation graph ==" not in s:
    needle = """echo
echo "== Inspect operation graph =="
./inspect_bigeye_v2_caa_operation_graph.sh"""

    replacement = """echo
echo "== Inspect operation graph =="
./inspect_bigeye_v2_caa_operation_graph.sh

echo
echo "== Validate operation graph =="
./validate_bigeye_v2_caa_operation_graph.sh"""

    if needle not in s:
        raise SystemExit("could not patch build_bigeye_v2_caa.sh")

    s = s.replace(needle, replacement, 1)

p.write_text(s)
print(f"updated {p}")
PY

echo "added CAA operation-graph validation"
echo
echo "Run:"
echo "  ./validate_bigeye_v2_caa_operation_graph.sh"
echo "  ./build_bigeye_v2_caa.sh"
echo "  ./run_bigeye_v2_regression_suite.sh"
