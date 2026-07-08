#!/usr/bin/env bash
set -euo pipefail

mkdir -p tools/caa

cat > tools/caa/README.md <<'MD'
# CAA Toolchain

The CAA Toolchain transforms declarative assessment descriptions into generated documentation, validation reports, and eventually runtime components.

The runtime does not depend on the toolchain.

## Inputs

- `caa_manifest.yml`
- `package.meta`

## Outputs

- package catalog
- package registry
- architecture summaries
- validation reports
- generated assessment-cycle code

## Principles

- Metadata is the single source of truth.
- Generated files are never edited by hand.
- Scientific code is independent of generators.
- Every generated artifact should be reproducible.

## Flow

```text
Manifest + Metadata
        |
        v
   CAA Toolchain
        |
        v
Generated Artifacts
        |
        v
 AssessmentCycle
        |
        v
 Quadra Solver
```
MD

cat > tools/caa/generate_catalog.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path

PACKAGE_ORDER = [
    "life_history",
    "population",
    "movement",
    "fleet",
    "observation",
    "likelihood",
]

BASE = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages")
OUT = Path("examples/NMFS/pifsc_bigeye_tuna/v2/CAA_PACKAGE_CATALOG.md")


def parse_meta(path: Path) -> dict:
    result = {"steps": []}
    in_steps = False

    for line in path.read_text().splitlines():
        line = line.rstrip()

        if line == "steps:":
            in_steps = True
            continue

        if in_steps and line.startswith("  - "):
            result["steps"].append(line.replace("  - ", "", 1))
            continue

        in_steps = False

        if ": " in line:
            key, value = line.split(": ", 1)
            result[key] = value

    return result


def main() -> None:
    lines = [
        "# CAA Package Catalog",
        "",
        "Generated from `architecture/packages/*/package.meta`.",
        "",
    ]

    for package_dir in PACKAGE_ORDER:
        meta_path = BASE / package_dir / "package.meta"
        if not meta_path.exists():
            raise SystemExit(f"missing package metadata: {meta_path}")

        meta = parse_meta(meta_path)

        lines.extend([
            f"## {meta['name']}",
            "",
            f"**Purpose:** {meta['purpose']}",
            "",
            f"**Consumes:** {meta['consumes']}",
            "",
            f"**Produces:** {meta['produces']}",
            "",
            "**Steps:**",
            "",
        ])

        for step in meta["steps"]:
            lines.append(f"- {step}")

        lines.extend(["", "---", ""])

    OUT.write_text("\n".join(lines))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
PY

chmod +x tools/caa/generate_catalog.py

cat > generate_bigeye_v2_caa_catalog.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

python3 tools/caa/generate_catalog.py
SH

chmod +x generate_bigeye_v2_caa_catalog.sh

echo "created CAA toolchain README and catalog generator"
echo
echo "Run:"
echo "  ./generate_bigeye_v2_caa_catalog.sh"
echo "  head -80 examples/NMFS/pifsc_bigeye_tuna/v2/CAA_PACKAGE_CATALOG.md"
