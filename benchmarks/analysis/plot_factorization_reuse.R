#!/usr/bin/env Rscript

csv_path <- "benchmarks/exact_laplace_gradient/factorization_reuse_benchmark.csv"

if (!file.exists(csv_path)) {
  stop("Factorization reuse benchmark CSV not found.")
}

df <- read.csv(csv_path)

summary <- aggregate(
  cbind(fresh_ms, reuse_ms, reuse_ratio) ~ n_state,
  data = df,
  FUN = mean
)

png(
  filename = "benchmarks/analysis/factorization_reuse_ratio.png",
  width = 1000,
  height = 700
)

plot(
  summary$n_state,
  summary$reuse_ratio,
  type = "b",
  pch = 1,
  lty = 1,
  xlab = "Number of latent states",
  ylab = "reuse_ms / fresh_ms",
  main = "Sparse Factorization Reuse Efficiency"
)

abline(h = 1.0, lty = 2)

dev.off()

png(
  filename = "benchmarks/analysis/factorization_reuse_times.png",
  width = 1000,
  height = 700
)

ylim <- range(c(summary$fresh_ms, summary$reuse_ms), na.rm = TRUE)

plot(
  summary$n_state,
  summary$fresh_ms,
  type = "b",
  pch = 1,
  lty = 1,
  ylim = ylim,
  xlab = "Number of latent states",
  ylab = "Factorization time (ms)",
  main = "Fresh vs Reused Sparse Factorization"
)

lines(
  summary$n_state,
  summary$reuse_ms,
  type = "b",
  pch = 2,
  lty = 2
)

legend(
  "topleft",
  legend = c("fresh", "reuse"),
  pch = c(1, 2),
  lty = c(1, 2)
)

dev.off()

png(
  filename = "benchmarks/analysis/factorization_reuse_vs_fill.png",
  width = 1000,
  height = 700
)

fill_summary <- aggregate(
  fill_ratio ~ n_state,
  data = df,
  FUN = mean
)

plot(
  fill_summary$fill_ratio,
  summary$reuse_ratio,
  type = "b",
  pch = 1,
  lty = 1,
  xlab = "Fill ratio",
  ylab = "reuse_ms / fresh_ms",
  main = "Reuse Efficiency vs Fill Ratio"
)

dev.off()

cat("Wrote factorization reuse plots.\n")
