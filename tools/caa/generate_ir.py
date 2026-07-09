#!/usr/bin/env python3
from pathlib import Path
import json

ORDER = ["life_history", "population", "movement", "fleet", "observation", "likelihood"]
BASE = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_IR.json")


def split_list(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def parse_meta(path: Path) -> dict:
    out = {"steps": [], "operations": []}
    in_steps = False
    in_operations = False

    for line in path.read_text().splitlines():
        line = line.rstrip()

        if line == "steps:":
            in_steps = True
            in_operations = False
            continue

        if line == "operations:":
            in_operations = True
            in_steps = False
            continue

        if in_steps and line.startswith("  - "):
            out["steps"].append(line[4:])
            continue

        if in_operations and line.startswith("  - "):
            out["operations"].append(line[4:])
            continue

        in_steps = False
        in_operations = False

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
        "operations": meta.get("operations", []),
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
