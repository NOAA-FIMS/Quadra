#!/usr/bin/env Rscript

csv_path <- "benchmarks/comparisons/tmb_state_space/comparison_outputs/quadra_state_space_compare.csv"

if (!file.exists(csv_path)) {
  stop("State-space comparison CSV not found.")
}

df <- read.csv(csv_path)

required <- c(
  "n_state",
  "hessian_nnz",
  "hessian_density",
  "factor_nnz",
  "fill_ratio",
  "workspace_ms",
  "factorization_ms"
)

missing <- setdiff(required, names(df))

if (length(missing) > 0) {
  stop(paste("Missing columns:", paste(missing, collapse = ", ")))
}

png(
  filename = "benchmarks/analysis/state_space_hessian_nnz.png",
  width = 1000,
  height = 700
)

plot(
  df$n_state,
  df$hessian_nnz,
  type = "b",
  xlab = "Number of latent states",
  ylab = "nnz(H_uu)",
  main = "State-Space Random-Effect Hessian Sparsity"
)

dev.off()

png(
  filename = "benchmarks/analysis/state_space_hessian_density.png",
  width = 1000,
  height = 700
)

plot(
  df$n_state,
  df$hessian_density,
  type = "b",
  xlab = "Number of latent states",
  ylab = "Hessian density",
  main = "State-Space Hessian Density"
)

dev.off()

png(
  filename = "benchmarks/analysis/state_space_factor_fill.png",
  width = 1000,
  height = 700
)

plot(
  df$n_state,
  df$fill_ratio,
  type = "b",
  xlab = "Number of latent states",
  ylab = "factor nnz / Hessian nnz",
  main = "State-Space Sparse Factor Fill Ratio"
)

dev.off()

png(
  filename = "benchmarks/analysis/state_space_timing_vs_structure.png",
  width = 1000,
  height = 700
)

plot(
  df$hessian_nnz,
  df$workspace_ms,
  type = "b",
  xlab = "nnz(H_uu)",
  ylab = "Workspace build time (ms)",
  main = "State-Space Timing vs Hessian Sparsity"
)

dev.off()

cat("Wrote state-space structure plots.\n")
