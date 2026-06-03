#!/usr/bin/env bash
set -euo pipefail

OUTDIR="benchmarks/state_space_surplus_scaling"

if [[ ! -d "$OUTDIR" ]]; then
  echo "ERROR: missing $OUTDIR"
  echo "Run install_quadra_tmb_scaling_artifact_v1.sh first."
  exit 1
fi

cat > "$OUTDIR/make_scaling_plot.R" <<'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("usage: make_scaling_plot.R results.csv output.png")
}

csv_path <- args[[1]]
out_path <- args[[2]]

d <- read.csv(csv_path)

png(out_path, width = 1200, height = 800, res = 160)

ylim <- range(c(d$quadra_ms, d$tmb_ms), finite = TRUE)

plot(
  d$n,
  d$quadra_ms,
  type = "b",
  log = "y",
  pch = 16,
  lwd = 2,
  ylim = ylim,
  xlab = "Number of years / latent-state scale",
  ylab = "Milliseconds per fixed-theta Laplace evaluation",
  main = "State-space surplus production Laplace scaling"
)

lines(d$n, d$tmb_ms, type = "b", pch = 17, lwd = 2)

legend(
  "topleft",
  legend = c("Quadra analytic tridiagonal", "TMB AD/Laplace"),
  pch = c(16, 17),
  lwd = 2,
  bty = "n"
)

grid()

dev.off()

cat("wrote", out_path, "\n")
EOF

chmod +x "$OUTDIR/make_scaling_plot.R"

target="run_state_space_surplus_scaling_artifact.sh"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/run_state_space_surplus_scaling_artifact.sh.plot_fallback.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("run_state_space_surplus_scaling_artifact.sh")
s = p.read_text()

old = '''if command -v python3 >/dev/null 2>&1; then
  python3 "$OUTDIR/make_scaling_plot.py" "$CSV" "$PLOT" || {
    echo "Plot generation failed; CSV is still available: $CSV"
  }
fi
'''

new = '''plot_done=0

if command -v python3 >/dev/null 2>&1; then
  if python3 "$OUTDIR/make_scaling_plot.py" "$CSV" "$PLOT"; then
    plot_done=1
  else
    echo "Python/matplotlib plot failed; trying R base graphics fallback..."
  fi
fi

if [[ "$plot_done" -eq 0 ]] && command -v Rscript >/dev/null 2>&1; then
  if Rscript "$OUTDIR/make_scaling_plot.R" "$CSV" "$PLOT"; then
    plot_done=1
  fi
fi

if [[ "$plot_done" -eq 0 ]]; then
  echo "Plot generation failed; CSV is still available: $CSV"
fi
'''

if old not in s:
    raise SystemExit("Could not find plotting block to patch")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Added R fallback plotting for the scaling artifact.

Run:
  ./run_state_space_surplus_scaling_artifact.sh 10 25,50,100,250

Or plot an existing CSV directly:
  Rscript benchmarks/state_space_surplus_scaling/make_scaling_plot.R \
    benchmarks/state_space_surplus_scaling/results.csv \
    benchmarks/state_space_surplus_scaling/scaling_plot.png

EOF
