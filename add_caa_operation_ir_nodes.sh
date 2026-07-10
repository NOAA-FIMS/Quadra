#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

python3 - <<'PY'
from pathlib import Path
import json

base = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")

nodes = {
    "life_history": [
        {"id": "ComputeNaturalMortality", "order": 0,
         "reads": ["LifeHistoryParameters.log_m_young_offset", "LifeHistoryParameters.log_m_old_offset"],
         "writes": ["LifeHistoryState.m_at_age"]},
        {"id": "ComputeWeightAtAge", "order": 1,
         "reads": ["BigeyeModelData"],
         "writes": ["LifeHistoryState.weight_at_age"]},
        {"id": "ComputeMaturityAtAge", "order": 2,
         "reads": ["BigeyeModelData"],
         "writes": ["LifeHistoryState.maturity_at_age"]},
    ],
    "population": [
        {"id": "Recruitment", "order": 0,
         "reads": ["PopulationParameters.r0"],
         "writes": ["PopulationState.recruits_by_year", "PopulationState.numbers_at_age"]},
        {"id": "Survival", "order": 1,
         "reads": ["PopulationState.numbers_at_age", "FleetState.z_at_age"],
         "writes": ["PopulationState.survivors_at_age"]},
        {"id": "Aging", "order": 2,
         "reads": ["PopulationState.survivors_at_age"],
         "writes": ["PopulationState.numbers_at_age"]},
        {"id": "PlusGroup", "order": 3,
         "reads": ["PopulationState.survivors_at_age", "PopulationState.numbers_at_age"],
         "writes": ["PopulationState.numbers_at_age"]},
        {"id": "SpawningBiomass", "order": 4,
         "reads": ["PopulationState.numbers_at_age", "LifeHistoryState.weight_at_age",
                   "LifeHistoryState.maturity_at_age"],
         "writes": ["PopulationState.spawning_biomass_by_year"]},
    ],
    "movement": [
        {"id": "IdentityMovement", "order": 0,
         "reads": ["PopulationState.numbers_at_age"],
         "writes": ["PopulationState.numbers_at_age"]},
    ],
    "fleet": [
        {"id": "LogisticSelectivity", "order": 0,
         "reads": ["FleetParameters.sel_a50", "FleetParameters.sel_slope"],
         "writes": ["FleetState.selectivity_at_age"]},
        {"id": "FishingMortality", "order": 1,
         "reads": ["FleetParameters.fbar", "FleetState.selectivity_at_age",
                   "LifeHistoryState.m_at_age"],
         "writes": ["FleetState.f_at_age", "FleetState.z_at_age"]},
        {"id": "BaranovCatch", "order": 2,
         "reads": ["PopulationState.numbers_at_age", "FleetState.f_at_age",
                   "FleetState.z_at_age", "LifeHistoryState.weight_at_age"],
         "writes": ["FleetState.catch_numbers_at_age", "FleetState.catch_biomass_at_age",
                    "FleetState.total_catch_biomass_by_year"]},
    ],
    "observation": [
        {"id": "BiomassIndexPrediction", "order": 0,
         "reads": ["PopulationState.spawning_biomass_by_year", "FleetParameters.q_index"],
         "writes": ["FleetState.predicted_index_by_year"]},
        {"id": "CatchAgeCompositionPrediction", "order": 1,
         "reads": ["FleetState.catch_numbers_at_age"],
         "writes": ["FleetState.predicted_catch_age_proportion"]},
    ],
    "likelihood": [
        {"id": "LognormalCatchLikelihood", "order": 0,
         "reads": ["BigeyeModelData.observed_catch_biomass_by_year",
                   "FleetState.total_catch_biomass_by_year", "FleetParameters.catch_sigma"],
         "writes": ["LikelihoodState.catch_nll"]},
        {"id": "LognormalIndexLikelihood", "order": 1,
         "reads": ["BigeyeModelData.observed_index_by_year",
                   "FleetState.predicted_index_by_year", "FleetParameters.index_sigma"],
         "writes": ["LikelihoodState.index_nll"]},
        {"id": "MultinomialAgeCompLikelihood", "order": 2,
         "reads": ["BigeyeModelData.observed_catch_age_proportion",
                   "BigeyeModelData.catch_agecomp_sample_size",
                   "FleetState.predicted_catch_age_proportion"],
         "writes": ["LikelihoodState.agecomp_nll"]},
        {"id": "SumLikelihood", "order": 3,
         "reads": ["LikelihoodState.catch_nll", "LikelihoodState.index_nll",
                   "LikelihoodState.agecomp_nll"],
         "writes": ["LikelihoodState.total_nll"]},
    ],
}

package_names = {
    "life_history": "LifeHistoryPackage",
    "population": "PopulationPackage",
    "movement": "MovementPackage",
    "fleet": "FleetPackage",
    "observation": "ObservationPackage",
    "likelihood": "LikelihoodPackage",
}

for key, operation_nodes in nodes.items():
    path = base / key / "operation_nodes.json"
    path.write_text(json.dumps({
        "package": package_names[key],
        "operations": operation_nodes,
    }, indent=2) + "\n")
    print(f"wrote {path}")
PY

python3 - <<'PY'
from pathlib import Path

p = Path("tools/caa/generate_ir.py")
s = p.read_text()

if "def load_operation_nodes" not in s:
    marker = "\n\npackages = []"
    helper = '''
def load_operation_nodes(meta_path: Path) -> list[dict]:
    node_path = meta_path.parent / "operation_nodes.json"
    if not node_path.exists():
        return []

    payload = json.loads(node_path.read_text())
    nodes = payload.get("operations", [])

    for node in nodes:
        node.setdefault("reads", [])
        node.setdefault("writes", [])
        node.setdefault("order", 0)

    return sorted(nodes, key=lambda node: (node["order"], node["id"]))
'''
    if marker not in s:
        raise SystemExit("could not find packages insertion point")
    s = s.replace(marker, "\n\n" + helper + marker, 1)

if '"operation_nodes": operation_nodes,' not in s:
    old = 'for key in ORDER:\n    meta = parse_meta(BASE / key / "package.meta")\n\n    package = {'
    new = 'for key in ORDER:\n    meta_path = BASE / key / "package.meta"\n    meta = parse_meta(meta_path)\n    operation_nodes = load_operation_nodes(meta_path)\n\n    package = {'
    if old not in s:
        raise SystemExit("could not patch package loop")
    s = s.replace(old, new, 1)

    old_field = '        "operations": meta.get("operations", []),'
    new_field = '        "operations": meta.get("operations", []),\n        "operation_nodes": operation_nodes,'
    if old_field not in s:
        raise SystemExit("could not patch operation_nodes field")
    s = s.replace(old_field, new_field, 1)

p.write_text(s)
print(f"updated {p}")
PY

cat > tools/caa/generate_operation_graph.py <<'PY'
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
PY

chmod +x tools/caa/generate_operation_graph.py

cat > generate_bigeye_v2_caa_operation_graph.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
./generate_bigeye_v2_caa_ir.sh
python3 tools/caa/generate_operation_graph.py
SH

chmod +x generate_bigeye_v2_caa_operation_graph.sh

cat > inspect_bigeye_v2_caa_operation_graph.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
from pathlib import Path
import json

graph = json.loads(Path(
    "examples/NMFS/pifsc_bigeye_tuna/v2/CAA_OPERATION_GRAPH.json"
).read_text())

print("CAA Operation Graph")
print()
print(f"nodes: {len(graph['nodes'])}")
print(f"edges: {len(graph['edges'])}")
print()

for node in graph["nodes"]:
    print(node["key"])
    print(f"  order: {node['order']}")
    print(f"  reads: {', '.join(node['reads']) or '(none)'}")
    print(f"  writes: {', '.join(node['writes']) or '(none)'}")
    print()

if graph["diagnostics"]:
    print("diagnostics:")
    for d in graph["diagnostics"]:
        print(f"  {d['level']}: {d['kind']} {d['operation']} <- {d['field']}")
else:
    print("diagnostics: clean")
PY
SH

chmod +x inspect_bigeye_v2_caa_operation_graph.sh

python3 - <<'PY'
from pathlib import Path

p = Path("build_bigeye_v2_caa.sh")
s = p.read_text()

if "== Generate operation graph ==" not in s:
    needle = 'echo\necho "== Generate IR =="\n./generate_bigeye_v2_caa_ir.sh'
    replacement = '''echo
echo "== Generate IR =="
./generate_bigeye_v2_caa_ir.sh

echo
echo "== Generate operation graph =="
./generate_bigeye_v2_caa_operation_graph.sh

echo
echo "== Inspect operation graph =="
./inspect_bigeye_v2_caa_operation_graph.sh'''
    if needle not in s:
        raise SystemExit("could not patch build_bigeye_v2_caa.sh")
    s = s.replace(needle, replacement, 1)

p.write_text(s)
print(f"updated {p}")
PY

echo
echo "Added CAA operation nodes and operation graph."
echo
echo "Run:"
echo "  ./generate_bigeye_v2_caa_operation_graph.sh"
echo "  ./inspect_bigeye_v2_caa_operation_graph.sh"
echo "  ./build_bigeye_v2_caa.sh"
echo "  ./run_bigeye_v2_regression_suite.sh"
