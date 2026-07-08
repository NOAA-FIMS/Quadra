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
    return out


packages = []
created = {}
updated = {}
consumed = {}

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
        "steps": meta["steps"],
    }

    packages.append(package)

    for item in package["creates"]:
        created.setdefault(item, []).append(package["name"])

    for item in package["updates"]:
        updated.setdefault(item, []).append(package["name"])

    for item in package["consumes"]:
        consumed.setdefault(item, []).append(package["name"])


diagnostics = []

for state, creators in created.items():
    if len(creators) > 1:
        diagnostics.append({
            "level": "warning",
            "kind": "multiple_creators",
            "state": state,
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

ir = {
    "name": "Bigeye v2 CAA IR",
    "packages": packages,
    "created": created,
    "updated": updated,
    "consumed": consumed,
    "diagnostics": diagnostics,
}

OUT.write_text(json.dumps(ir, indent=2) + "\n")
print(f"wrote {OUT}")

if diagnostics:
    print("diagnostics:")
    for d in diagnostics:
        print(f"  {d['level']}: {d['kind']} {d.get('state', '')}")
else:
    print("diagnostics: clean")
