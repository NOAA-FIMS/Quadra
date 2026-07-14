#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
from pathlib import Path

p = Path("tools/caa/validate_operation_graph.py")
s = p.read_text()

# Replace the old topological cycle block with planner-aware SCC handling.
start_marker = "# ---------------------------------------------------------------------------\n# Cycle detection\n# ---------------------------------------------------------------------------"
end_marker = "# ---------------------------------------------------------------------------\n# Declared-order contradiction checks\n# ---------------------------------------------------------------------------"

start = s.find(start_marker)
end = s.find(end_marker)
if start == -1 or end == -1 or end <= start:
    raise SystemExit("could not locate cycle-validation block")

replacement = '''# ---------------------------------------------------------------------------
# Planner-aware feedback-cycle validation
#
# The static operation graph can contain legitimate feedback. In Bigeye,
# PopulationState feeds FleetPackage and FleetState feeds PopulationPackage.
# The execution planner resolves this with a FleetPackage bootstrap pass and
# a later rerun. Such SCCs are warnings when they contain a planner-declared
# rerunnable package; all other SCCs remain errors.
# ---------------------------------------------------------------------------
PLAN = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_EXECUTION_PLAN.json")
plan = json.loads(PLAN.read_text()) if PLAN.exists() else {}
rerunnable_packages = set(plan.get("rerunnable_packages", []))

index_counter = [0]
stack = []
on_stack = set()
indices = {}
lowlink = {}
components = []


def strongconnect(node):
    indices[node] = index_counter[0]
    lowlink[node] = index_counter[0]
    index_counter[0] += 1

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

    packages = sorted({
        node_by_key[key]["package"]
        for key in component
        if key in node_by_key
    })
    rerunnable_in_component = sorted(set(packages) & rerunnable_packages)

    fields = set()
    component_set = set(component)
    for edge in edges:
        if (
            edge.get("from") in component_set
            and edge.get("to") in component_set
        ):
            fields.update(edge.get("fields", []))

    details = {
        "operations": component,
        "packages": packages,
        "fields": sorted(fields),
        "rerunnable_packages": rerunnable_in_component,
        "resolution": (
            "execution planner bootstrap/rerun"
            if rerunnable_in_component
            else "unresolved"
        ),
    }

    if rerunnable_in_component:
        warnings.append({
            "level": "warning",
            "kind": "operation_graph_feedback_cycle",
            **details,
        })
    else:
        errors.append({
            "level": "error",
            "kind": "operation_graph_cycle",
            **details,
        })


'''

s = s[:start] + replacement + s[end:]
p.write_text(s)
print(f"updated {p}")
PY

python3 - <<'PY'
from pathlib import Path

p = Path("validate_bigeye_v2_caa_operation_graph.sh")
s = p.read_text()

old = '''./generate_bigeye_v2_caa_operation_graph.sh
python3 tools/caa/validate_operation_graph.py'''
new = '''./generate_bigeye_v2_caa_operation_graph.sh
./generate_bigeye_v2_caa_execution_plan.sh
python3 tools/caa/validate_operation_graph.py'''

if old in s:
    s = s.replace(old, new, 1)
elif "./generate_bigeye_v2_caa_execution_plan.sh" not in s:
    raise SystemExit("could not patch operation graph validation wrapper")

p.write_text(s)
print(f"updated {p}")
PY

chmod +x tools/caa/validate_operation_graph.py
chmod +x validate_bigeye_v2_caa_operation_graph.sh

echo "refined operation-graph cycle validation"
echo
echo "Run:"
echo "  ./validate_bigeye_v2_caa_operation_graph.sh"
echo "  ./build_bigeye_v2_caa.sh"
echo "  ./run_bigeye_v2_regression_suite.sh"
