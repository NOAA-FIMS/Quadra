#!/usr/bin/env python3
"""Validate Quadra JSONL benchmark artifacts without third-party packages."""

import json
import math
import pathlib
import sys

REQUIRED = {
    "schema_version", "benchmark", "model_form", "phase", "sample",
    "dimension", "hessian_nnz", "backend", "elapsed_ms", "success",
}
MODEL_FORMS = {
    "diagonal", "tridiagonal", "banded", "block_diagonal", "arrowhead",
    "general_sparse", "nearly_dense", "dense", "parameter_dependent",
}
PHASES = {
    "record", "replay", "gradient", "hessian", "mode_solve", "pattern_discovery",
    "symbolic_analysis", "factorization", "logdet", "hdot", "cold_total",
    "warm_total",
}


def fail(path, line_number, message):
    raise ValueError(f"{path}:{line_number}: {message}")


def validate(path):
    count = 0
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            count += 1
            try:
                row = json.loads(line)
            except json.JSONDecodeError as error:
                fail(path, line_number, f"invalid JSON: {error}")
            missing = REQUIRED - row.keys()
            if missing:
                fail(path, line_number, f"missing fields: {sorted(missing)}")
            if row["schema_version"] != 1:
                fail(path, line_number, "unsupported schema_version")
            if row["model_form"] not in MODEL_FORMS:
                fail(path, line_number, "unknown model_form")
            if row["phase"] not in PHASES:
                fail(path, line_number, "unknown phase")
            for field in ("sample", "dimension", "hessian_nnz"):
                if not isinstance(row[field], int) or row[field] < 0:
                    fail(path, line_number, f"{field} must be nonnegative integer")
            elapsed = row["elapsed_ms"]
            if not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed < 0:
                fail(path, line_number, "elapsed_ms must be finite and nonnegative")
            if not isinstance(row["success"], bool):
                fail(path, line_number, "success must be boolean")
    if count == 0:
        raise ValueError(f"{path}: artifact contains no results")
    return count


def main(argv):
    if len(argv) < 2:
        raise SystemExit("usage: validate_results.py RESULTS.jsonl [...]")
    total = 0
    for name in argv[1:]:
        total += validate(pathlib.Path(name))
    print(f"validated {total} benchmark results")


if __name__ == "__main__":
    main(sys.argv)
