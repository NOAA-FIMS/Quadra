#!/usr/bin/env python3

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NORMALIZED = ROOT / "normalized"
OUTPUTS = ROOT / "outputs"

def read_csv(path):
    if not path.exists():
        return []

    with open(path, newline="") as f:
        return list(csv.DictReader(f))

def to_float(x):
    try:
        if x is None or x == "":
            return None
        return float(x)
    except ValueError:
        return None

def summarize_engine(rows, x_col, y_col):
    clean = []
    for r in rows:
        x = to_float(r.get(x_col))
        y = to_float(r.get(y_col))
        if x is not None and y is not None:
            clean.append((x, y))

    if len(clean) < 2:
        return None

    clean.sort()
    x0, y0 = clean[0]
    x1, y1 = clean[-1]

    ratio = y1 / y0 if y0 not in (None, 0.0) else None

    return {
        "x0": x0,
        "y0": y0,
        "x1": x1,
        "y1": y1,
        "ratio": ratio,
    }

def fmt(x):
    if x is None:
        return "NA"
    if abs(x) >= 100:
        return f"{x:.1f}"
    if abs(x) >= 1:
        return f"{x:.3f}"
    return f"{x:.6f}"

def write_section_random_intercept(out, rows):
    out.append("## Random Intercept Benchmark")
    out.append("")

    if not rows:
        out.append("_No random-intercept benchmark rows found._")
        out.append("")
        return

    engines = sorted(set(r.get("engine", "") for r in rows if r.get("engine")))

    for engine in engines:
        engine_rows = [r for r in rows if r.get("engine") == engine]

        summary = summarize_engine(
            engine_rows,
            "n_obs",
            "total_wall_ms")

        if summary is None:
            out.append(f"- **{engine}**: insufficient total wall-time data.")
            continue

        out.append(
            f"- **{engine}** total wall time increased from "
            f"{fmt(summary['y0'])} ms at n={fmt(summary['x0'])} "
            f"to {fmt(summary['y1'])} ms at n={fmt(summary['x1'])} "
            f"({fmt(summary['ratio'])}x)."
        )

    out.append("")

    out.append("Interpretation notes:")
    out.append("")
    out.append("- Random-intercept benchmarks are useful for checking overhead and basic scaling.")
    out.append("- They are not the strongest test of sparse structure because there is only one random effect.")
    out.append("- Decomposed timings should be interpreted by phase, not as a single headline comparison.")
    out.append("")

def write_section_state_space(out, rows):
    out.append("## State-Space Benchmark")
    out.append("")

    if not rows:
        out.append("_No state-space benchmark rows found._")
        out.append("")
        return

    engines = sorted(set(r.get("engine", "") for r in rows if r.get("engine")))

    for engine in engines:
        engine_rows = [r for r in rows if r.get("engine") == engine]

        time_summary = summarize_engine(
            engine_rows,
            "n_state",
            "total_wall_ms")

        if time_summary is not None:
            out.append(
                f"- **{engine}** total wall time increased from "
                f"{fmt(time_summary['y0'])} ms at n_state={fmt(time_summary['x0'])} "
                f"to {fmt(time_summary['y1'])} ms at n_state={fmt(time_summary['x1'])} "
                f"({fmt(time_summary['ratio'])}x)."
            )

        nnz_summary = summarize_engine(
            engine_rows,
            "n_state",
            "hessian_nnz")

        if nnz_summary is not None:
            out.append(
                f"- **{engine}** Hessian nonzeros increased from "
                f"{fmt(nnz_summary['y0'])} to {fmt(nnz_summary['y1'])} "
                f"over the same state-size range."
            )

        fill_values = [
            to_float(r.get("fill_ratio"))
            for r in engine_rows
            if to_float(r.get("fill_ratio")) is not None
        ]

        if fill_values:
            out.append(
                f"- **{engine}** fill ratio ranged from "
                f"{fmt(min(fill_values))} to {fmt(max(fill_values))}."
            )

    out.append("")
    out.append("Interpretation notes:")
    out.append("")
    out.append("- State-space benchmarks are more informative for sparse mixed-effects scaling.")
    out.append("- Hessian nonzeros, factor nonzeros, and fill ratio help explain timing changes.")
    out.append("- Runtime should be interpreted alongside structure metrics, not separately.")
    out.append("")

def main():
    random_rows = read_csv(NORMALIZED / "random_intercept_normalized.csv")
    state_rows = read_csv(NORMALIZED / "state_space_normalized.csv")

    out = []
    out.append("# Quadra Benchmark Interpretation Summary")
    out.append("")
    out.append("This summary is generated automatically from normalized benchmark CSV outputs.")
    out.append("")

    write_section_random_intercept(out, random_rows)
    write_section_state_space(out, state_rows)

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    path = OUTPUTS / "benchmark_summary.md"
    path.write_text("\n".join(out) + "\n")

    print(f"Wrote benchmark summary: {path}")

if __name__ == "__main__":
    main()
