#!/usr/bin/env python3
"""Run each Laplace model form in a fresh process and combine its JSONL."""

import argparse
import pathlib
import subprocess
import sys

MODEL_FORMS = (
    "diagonal",
    "tridiagonal",
    "banded",
    "block_diagonal",
    "general_sparse",
    "dense",
)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", default=(
        "benchmarks/framework/laplace_model_catalog_benchmark"
    ))
    parser.add_argument("--dimension", type=int, default=128)
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv[1:])

    executable = pathlib.Path(args.executable)
    if not executable.is_file():
        raise SystemExit(f"benchmark executable not found: {executable}")
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as stream:
        for model_form in MODEL_FORMS:
            completed = subprocess.run(
                [str(executable), str(args.dimension), str(args.repetitions),
                 model_form],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
            )
            stream.write(completed.stdout)

    print(f"wrote {output} from {len(MODEL_FORMS)} fresh processes")


if __name__ == "__main__":
    main(sys.argv)
