#!/usr/bin/env bash
set -euo pipefail

# install_quadra_tmb_scaling_artifact_v1.sh
#
# Creates a clean benchmark artifact for the state-space surplus production
# scaling result:
#
#   - benchmark README
#   - CSV capture script
#   - Python plotting script
#   - one-shot run script
#
# This is meant to turn the current terminal benchmark into something
# reproducible and shareable.

mkdir -p benchmarks/state_space_surplus_scaling

cat > benchmarks/state_space_surplus_scaling/README.md <<'EOF'
# State-space surplus production scaling benchmark

This benchmark compares fixed-theta Laplace objective evaluation for the same
state-space surplus production model in:

```text
Quadra analytic latent-state tridiagonal implementation
TMB AD/Laplace implementation
```

The model is:

```text
pred_B[t+1] = B[t] + r B[t] (1 - B[t] / K) - C[t]
log_B[t+1]  = log(pred_B[t+1]) + process error
log(I[t])   = log(q) + log_B[t] + observation error
```

The fixed effects are held constant:

```text
r              = 0.5
K              = 700
q              = 0.0024
sigma_process  = 0.15
sigma_index    = 0.10
B0/K           = 0.90
```

The benchmark evaluates the marginal negative log likelihood after integrating
over latent log-biomass states.

## Current result

Representative run:

| n | Quadra ms/eval | TMB ms/eval | Quadra speedup |
|---:|---:|---:|---:|
| 25 | 0.027 | 0.100 | 3.7x |
| 50 | 0.060 | 0.600 | 9.9x |
| 100 | 0.117 | 2.200 | 18.7x |
| 250 | 0.255 | 14.100 | 55.2x |

Objectives matched to numerical precision.

## Interpretation

This benchmark does **not** show that Quadra is universally faster than TMB.

It shows that when the model has known latent Markov structure and Quadra
exploits the analytic tridiagonal Hessian, the Laplace evaluation can scale
nearly linearly and substantially outperform a generic AD/Laplace path.

## Run

From the repository root:

```bash
./run_state_space_surplus_scaling_artifact.sh 10 25,50,100,250
```

Outputs are written to:

```text
benchmarks/state_space_surplus_scaling/results.csv
benchmarks/state_space_surplus_scaling/scaling_plot.png
```
EOF

cat > benchmarks/state_space_surplus_scaling/make_scaling_plot.py <<'EOF'
#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> int:
    if len(sys.argv) < 3:
        print("usage: make_scaling_plot.py results.csv output.png")
        return 1

    csv_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    rows = []
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    quadra = []
    tmb = []

    for row in rows:
        n = int(row["n"])
        q = float(row["quadra_ms"])
        t = float(row["tmb_ms"])
        quadra.append((n, q))
        tmb.append((n, t))

    quadra.sort()
    tmb.sort()

    plt.figure()
    plt.plot([x for x, _ in quadra], [y for _, y in quadra], marker="o", label="Quadra analytic tridiagonal")
    plt.plot([x for x, _ in tmb], [y for _, y in tmb], marker="o", label="TMB AD/Laplace")
    plt.xlabel("Number of years / latent states scale")
    plt.ylabel("Milliseconds per fixed-theta Laplace evaluation")
    plt.yscale("log")
    plt.title("State-space surplus production Laplace scaling")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)

    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
EOF

chmod +x benchmarks/state_space_surplus_scaling/make_scaling_plot.py

cat > run_state_space_surplus_scaling_artifact.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250}"

OUTDIR="benchmarks/state_space_surplus_scaling"
mkdir -p "$OUTDIR"

QUADRA_LOG="$OUTDIR/quadra_raw.txt"
TMB_LOG="$OUTDIR/tmb_raw.txt"
CSV="$OUTDIR/results.csv"
PLOT="$OUTDIR/scaling_plot.png"

echo "Running Quadra benchmark..."
./run_quadra_scaled_analytic_latent_tridiagonal_benchmark.sh "$REPS" "$LENGTHS" | tee "$QUADRA_LOG"

echo
echo "Running TMB benchmark..."
./run_tmb_scaled_fixed_theta_benchmark.sh "$REPS" "$LENGTHS" | tee "$TMB_LOG"

python3 - <<'PYEOF'
from pathlib import Path
import csv
import re

outdir = Path("benchmarks/state_space_surplus_scaling")
quadra_log = outdir / "quadra_raw.txt"
tmb_log = outdir / "tmb_raw.txt"
csv_path = outdir / "results.csv"

quadra = {}
tmb = {}

q_line = re.compile(
    r"^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$"
)

t_line = re.compile(
    r"^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$"
)

for line in quadra_log.read_text().splitlines():
    m = q_line.match(line)
    if not m:
        continue
    n = int(m.group(1))
    quadra[n] = {
        "objective": float(m.group(2)),
        "joint": float(m.group(3)),
        "logdet": float(m.group(4)),
        "nnz": int(m.group(5)),
        "grad_norm": float(m.group(6)),
        "quadra_ms": float(m.group(7)),
    }

for line in tmb_log.read_text().splitlines():
    m = t_line.match(line)
    if not m:
        continue
    n = int(m.group(1))
    tmb[n] = {
        "tmb_objective": float(m.group(2)),
        "tmb_ms": float(m.group(3)),
    }

ns = sorted(set(quadra) & set(tmb))

with csv_path.open("w", newline="") as f:
    fieldnames = [
        "n",
        "quadra_objective",
        "tmb_objective",
        "objective_diff",
        "quadra_ms",
        "tmb_ms",
        "quadra_speedup",
        "joint",
        "logdet",
        "nnz",
        "grad_norm",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    for n in ns:
        q = quadra[n]
        t = tmb[n]
        qms = q["quadra_ms"]
        tms = t["tmb_ms"]
        writer.writerow({
            "n": n,
            "quadra_objective": q["objective"],
            "tmb_objective": t["tmb_objective"],
            "objective_diff": q["objective"] - t["tmb_objective"],
            "quadra_ms": qms,
            "tmb_ms": tms,
            "quadra_speedup": tms / qms if qms > 0 else "",
            "joint": q["joint"],
            "logdet": q["logdet"],
            "nnz": q["nnz"],
            "grad_norm": q["grad_norm"],
        })

print(f"wrote {csv_path}")
PYEOF

if command -v python3 >/dev/null 2>&1; then
  python3 "$OUTDIR/make_scaling_plot.py" "$CSV" "$PLOT" || {
    echo "Plot generation failed; CSV is still available: $CSV"
  }
fi

echo
echo "Benchmark artifact complete."
echo "CSV:  $CSV"
echo "Plot: $PLOT"
EOF

chmod +x run_state_space_surplus_scaling_artifact.sh

cat <<'EOF'

Installed state-space surplus scaling artifact.

Run:
  ./run_state_space_surplus_scaling_artifact.sh 10 25,50,100,250

Outputs:
  benchmarks/state_space_surplus_scaling/results.csv
  benchmarks/state_space_surplus_scaling/scaling_plot.png
  benchmarks/state_space_surplus_scaling/README.md

EOF
